import re
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolAlign, rdFMCS, Draw
import pandas as pd
import numpy as np
from rdkit.Chem.Draw import rdMolDraw2D

# ─────────────────────────────────────────
#  1) konwersja PDBQT → czysty PDB (bez pustych linii)
# ─────────────────────────────────────────
def convert_pdbqt_to_clean_pdb(src: Path, dst: Path):
    """
    Zamienia .pdbqt → .pdb:
    • bierze tylko linie ATOM/HETATM
    • rekalkuluje symbol pierwiastka (kolumny 77‑78)
    • NIE zostawia podwójnych \n, dodaje 'END' na końcu
    """
    out_lines = []
    with open(src) as fh:
        for ln in fh:
            if not ln.startswith(("ATOM", "HETATM")):
                continue
            if len(ln) < 78:
                continue

            atom_name = ln[12:16].strip()
            element = re.sub(r"[^A-Za-z]", "", atom_name).upper()  # np. 'C', 'CL'
            ln = ln[:76] + f"{element:>2}" + ln[78:]
            out_lines.append(ln.rstrip("\n"))  # usuwamy \n z oryginału

    if not out_lines:
        raise ValueError(f"No ATOM/HETATM lines found in {src}")

    out_lines.append("END")
    dst.write_text("\n".join(out_lines) + "\n")


# ─────────────────────────────────────────
#  2) „inteligentne” RMSD
# ─────────────────────────────────────────
def compute_rmsd(ref_pdb: Path, pred_pdb: Path) -> float:
    """
    RMSD heavy‑atomów z symetrią i dopasowaniem mapa‑do‑mapy:
    • usuwa hydrogeny
    • znajduje MCS i buduje atomMap
    • liczy 'best RMS' (uwzględnia permutacje symetryczne)
    """
    ref  = Chem.MolFromPDBFile(str(ref_pdb), removeHs=False, sanitize=False)
    pred = Chem.MolFromPDBFile(str(pred_pdb), removeHs=False, sanitize=False)
    if ref is None or pred is None:
        raise ValueError("RDKit could not parse molecules")

    ref_hvy  = Chem.RemoveHs(ref)
    pred_hvy = Chem.RemoveHs(pred)

    # 2a) różna liczba atomów → użyj MCS do mapowania
    if ref_hvy.GetNumAtoms() != pred_hvy.GetNumAtoms():
        mcs = rdFMCS.FindMCS(
            [ref_hvy, pred_hvy],
            completeRingsOnly=True,
            ringMatchesRingOnly=True,
            matchValences=True,
        )
        patt = Chem.MolFromSmarts(mcs.smartsString)
        ref_match  = ref_hvy.GetSubstructMatch(patt)
        pred_match = pred_hvy.GetSubstructMatch(patt)
        atom_map = list(zip(pred_match, ref_match))
        rmsd = rdMolAlign.AlignMol(pred_hvy, ref_hvy, atomMap=atom_map)
    else:
        # 2b) identyczna liczba atomów → spróbuj „najlepszego” RMSD
        rmsd = rdMolAlign.GetBestRMS(pred_hvy, ref_hvy)

    return rmsd



def draw_public_overlay(ref, pred, out_path, w=900, h=900):

    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    opts   = drawer.drawOptions()

    opts.padding         = 0.1
    opts.fixedBondLength = 22
    opts.bondLineWidth   = 1.0          # kontur natywnego

    # ── bezpieczne ustawienia czcionek ─────────────────────────────
    if hasattr(opts, "fontSize"):
        opts.fontSize = 14              # globalny rozmiar labeli
    if hasattr(opts, "minFontSize"):
        opts.minFontSize = 10

    opts.addAtomIndices = False
    opts.useBWAtomPalette()              # jednolite szarości

    # ── natywny (szary) ───────────────────────────────────────────
    grey = (0.3, 0.3, 0.3, 0.2)
    hl_ref = {i: grey for i in range(ref.GetNumAtoms())}
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer, ref,
        highlightAtoms=list(hl_ref),
        highlightAtomColors=hl_ref,
        highlightBondColors=hl_ref,
    )

    # ── dokowany (czerwony, grubszy, alfa) ─────────────────────────
    opts.bondLineWidth = 2.0
    red = (0.90, 0.15, 0.15, 0.6)         # RGBA
    hl_pred = {i: red for i in range(pred.GetNumAtoms())}
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer, pred,
        highlightAtoms=list(hl_pred),
        highlightAtomColors=hl_pred,
        highlightBondColors=hl_pred,
    )

    # ── zapis ──────────────────────────────────────────────────────
    drawer.FinishDrawing()
    Path(out_path).write_text(drawer.GetDrawingText())

def run_rmsd_and_plot(cfg, log):
    lig_dir = Path(cfg["paths"]["native_ligands_folder"])
    out_dir = Path(cfg["paths"]["output_folder"])
    results = []

    for docked_file in out_dir.glob("*__*__native_redock.pdbqt"):
        try:
            rec_id, lig_id, *_ = docked_file.stem.split("__")
            stem = f"{rec_id}_{lig_id}"
            native_pdb = lig_dir / f"{stem}.pdb"
            if not native_pdb.exists():
                log.warning(f"Native PDB not found for {stem}")
                continue

            converted_pdb = docked_file.with_suffix(".converted_fixed.pdb")
            if not converted_pdb.exists():
                convert_pdbqt_to_clean_pdb(docked_file, converted_pdb)

            # === RMSD (heavy atoms, best‑permutation) ===
            rmsd = compute_rmsd(native_pdb, converted_pdb)
            results.append((stem, rmsd))
            log.info(f"RMSD for {stem}: {rmsd:.2f} Å")

            # === Wizualizacje ===
            ref_vis  = Chem.RemoveHs(Chem.MolFromPDBFile(str(native_pdb)))
            pred_vis = Chem.RemoveHs(Chem.MolFromPDBFile(str(converted_pdb)))

            # 1. stary grid (zostawiamy, bo bywa przydatny)
            Draw.MolsToGridImage(
                [ref_vis, pred_vis],
                legends=["Native", "Docked"],
                molsPerRow=2,
                subImgSize=(300, 300),
            ).save(out_dir / f"{stem}_superposition.png")

            # 2. overlay
            # tworzymy grafikę (bez wodoru, żeby było czytelniej)
            ref_vis = Chem.RemoveHs(Chem.MolFromPDBFile(str(native_pdb)))
            pred_vis = Chem.RemoveHs(Chem.MolFromPDBFile(str(converted_pdb)))


            draw_public_overlay(ref_vis, pred_vis, out_dir / f"{stem}_overlay.svg")

        except Exception as e:
            log.error(f"Failed RMSD for {docked_file.name}: {e}")

    pd.DataFrame(results, columns=["Ligand", "RMSD"]).to_csv(
        out_dir / "rmsd_summary.csv", index=False
    )



