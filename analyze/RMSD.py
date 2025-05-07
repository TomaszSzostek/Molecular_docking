from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolAlign, Draw
import matplotlib.pyplot as plt
import pandas as pd

def compute_rmsd(ref_path: Path, pred_path: Path) -> float:
    ref = Chem.MolFromPDBFile(str(ref_path), removeHs=False)
    with open(pred_path, "r") as f:
        pdbqt_block = f.read()
    pred = Chem.MolFromPDBBlock(pdbqt_block, removeHs=False)

    if ref is None or pred is None:
        raise ValueError(f"Could not parse molecules: {ref_path} / {pred_path}")

    rmsd = rdMolAlign.AlignMol(pred, ref)
    return rmsd

def run_rmsd_and_plot(cfg, log):
    lig_dir = Path(cfg["paths"]["native_ligands_folder"])
    out_dir = Path(cfg["paths"]["output_folder"])
    results = []

    # Obsługa nazw w formacie {recID}__{ligID}__native_redock.pdbqt
    for docked_file in out_dir.glob("*__*__native_redock.pdbqt"):
        try:
            parts = docked_file.stem.split("__")
            if len(parts) < 3:
                continue

            rec_id, lig_id, mode = parts[0], parts[1], parts[2]
            stem = f"{rec_id}_{lig_id}"
            native_pdb = lig_dir / f"{stem}.pdb"

            if not native_pdb.exists():
                log.warning(f"Native PDB not found for {stem}")
                continue

            rmsd = compute_rmsd(native_pdb, docked_file)
            results.append((stem, rmsd))
            log.info(f"RMSD for {stem}: {rmsd:.2f}")

            # Wizualizacja 2D
            ref = Chem.MolFromPDBFile(str(native_pdb), removeHs=False)
            # Fallback: parsujemy tylko linie ATOM/HETATM i zapisujemy jako tymczasowy .pdb
            pdb_lines = [l for l in docked_file.read_text().splitlines() if
                         l.startswith("ATOM") or l.startswith("HETATM")]
            tmp_pdb = out_dir / f"{stem}_docked.pdb"
            tmp_pdb.write_text("\n".join(pdb_lines) + "\n")

            pred = Chem.MolFromPDBFile(str(tmp_pdb), removeHs=False)

            img = Draw.MolsToGridImage([ref, pred], legends=["Native", "Docked"])
            img.save(out_dir / f"{stem}_superposition.png")

        except Exception as e:
            log.error(f"Failed RMSD for {docked_file.name}: {e}")

    df = pd.DataFrame(results, columns=["Ligand", "RMSD"])
    df.to_csv(out_dir / "rmsd_summary.csv", index=False)



