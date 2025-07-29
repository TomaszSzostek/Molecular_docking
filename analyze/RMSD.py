from __future__ import annotations
"""
RMSD.py

Compute heavy‑atom RMSD and generate 2D overlay images between docked and native ligands.

This module implements:
  1. Direct Kabsch RMSD calculation on PDBQT coordinates to avoid ring/bond-order issues.
  2. Optional outlier refinement when RMSD exceeds threshold.
  3. Conversion of PDBQT to clean PDB for RDKit drawing.
  4. RDKit-based 2D overlay generation saved as SVG/PNG/PDF.
  5. Summary CSV with columns: RMSD_refined, n_common, n_common_refined.

Public API:
    run_rmsd_and_plot(cfg: dict, log: logging.Logger)

Configuration:
    cfg["paths"]["native_ligands_folder"] : Path to native ligands (.pdb or .pdbqt)
    cfg["paths"]["output_folder"]         : Path where docking PDBQT outputs live
    cfg["paths"]["visuals"]               : Path to write overlays and summary CSV
"""

from pathlib import Path
import re
import logging
from typing import List, Tuple, Dict

import numpy as np
import pandas as pd
import cairosvg

# ─────────────────────── RDKit only for drawing ───────────────────────
try:
    from rdkit import RDLogger
    if RDLogger is not None:
        RDLogger.DisableLog("rdApp.warning")  # suppress RDKit warnings
except Exception:
    RDLogger = None

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

# ───────────────────────────── Constants ─────────────────────────────────
MIN_COMMON_ATOMS = 4          # Minimum number of shared atoms for RMSD
REFINE_KEEP_PERCENT = 80      # Percent of atom pairs to keep after outlier removal

_PERIODIC = {
    "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
    "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
    "Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy",
    "Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb",
    "Bi","Po","At","Rn",
}

# ─────────────────────── Coordinates & RMSD (no RDKit) ────────────────────

def _guess_element(atom_name: str) -> str:
    """
    Infer chemical element from a PDB atom name (columns 13–16).

    Strips leading digits, handles two-letter elements only if second char is lowercase.
    """
    name = atom_name.strip()
    if not name:
        return "X"
    if name[0].isdigit():
        name = name[1:]  # remove leading digit
    # two-letter element if second char lowercase and valid
    if len(name) >= 2 and name[1].islower():
        two = name[:2].capitalize()
        if two in _PERIODIC:
            return two
    # fallback to first letter
    return name[0].upper()


def _read_coords(path: Path, heavy_only: bool = True):
    """
    Read atom names, elements, and XYZ coordinates from a PDB/PDBQT file.

    Parameters
    ----------
    path : Path
        Input file containing ATOM/HETATM lines.
    heavy_only : bool
        If True, skip hydrogen atoms.

    Returns
    -------
    (names, elems, xyz)
      names : np.ndarray of atom name strings
      elems : np.ndarray of element symbols
      xyz   : np.ndarray of shape (N,3) with float coordinates

    Raises
    ------
    ValueError if no atoms found.
    """
    names: list[str] = []
    elems: list[str] = []
    xyz:   list[Tuple[float,float,float]] = []
    for ln in path.read_text().splitlines():
        if ln.startswith(("ATOM","HETATM")):
            an = ln[12:16]
            el = _guess_element(an)
            if heavy_only and el == "H":
                continue
            try:
                x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            except ValueError:
                continue
            names.append(an.strip())
            elems.append(el)
            xyz.append((x,y,z))
    if not xyz:
        raise ValueError(f"No ATOM/HETATM in {path.name}")
    return np.array(names), np.array(elems), np.asarray(xyz, float)


def _map_by_name(ref_names: np.ndarray, pred_names: np.ndarray):
    """
    Map atoms by name between reference and predicted sets.

    Each occurrence is paired FIFO.
    """
    from collections import defaultdict, deque
    idxs = defaultdict(deque)
    for i,nm in enumerate(ref_names):
        idxs[nm].append(i)
    amap: list[Tuple[int,int]] = []
    for i,nm in enumerate(pred_names):
        if idxs[nm]:
            amap.append((i, idxs[nm].popleft()))
    return amap


def _kabsch(P: np.ndarray, Q: np.ndarray):
    """
    Perform Kabsch alignment between coordinate sets P and Q.

    Returns (rmsd, distances_per_atom).
    """
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    C  = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(C)
    if np.linalg.det(V)*np.linalg.det(Wt) < 0:
        V[:,-1] *= -1
    U   = V @ Wt
    Prot = Pc @ U
    dists = np.linalg.norm(Prot - Qc, axis=1)
    rmsd  = float(np.sqrt((dists**2).sum() / P.shape[0]))
    return rmsd, dists


def _refine_map_by_distance(ref_xyz: np.ndarray, pred_xyz: np.ndarray,
                            amap: List[Tuple[int,int]],
                            keep_percent=REFINE_KEEP_PERCENT):
    """
    Remove outlier atom pairs based on distance percentile threshold.
    """
    p_idx, r_idx = zip(*amap)
    P = pred_xyz[list(p_idx)]
    Q = ref_xyz[list(r_idx)]
    _, d = _kabsch(P, Q)
    thr = np.percentile(d, keep_percent)
    return [pair for k, pair in enumerate(amap) if d[k] <= thr]


def calc_rmsd_from_pdbqt(ref_fp: Path, pred_fp: Path,
                         heavy_only: bool = True) -> Dict[str, float|int|None]:
    """
    Compute RMSD metrics between native and predicted poses.

    Returns dict with keys:
      rmsd_full, n_common, n_ref, rmsd_refined, n_common_refined
    """
    ref_names, _, ref_xyz     = _read_coords(ref_fp, heavy_only)
    pred_names, _, pred_xyz   = _read_coords(pred_fp, heavy_only)

    amap = _map_by_name(ref_names, pred_names)
    if len(amap) < MIN_COMMON_ATOMS:
        if len(ref_xyz) == len(pred_xyz):
            amap = [(i,i) for i in range(len(ref_xyz))]
        else:
            raise ValueError(f"Too few common atoms ({len(amap)}) between {ref_fp.name} and {pred_fp.name}")

    # full RMSD
    p_idx, r_idx = zip(*amap)
    rmsd_full, _ = _kabsch(pred_xyz[list(p_idx)], ref_xyz[list(r_idx)])

    # optional refined RMSD
    rmsd_refined = None
    n_common_refined = None
    if rmsd_full > 2.0:
        amap_ref = _refine_map_by_distance(ref_xyz, pred_xyz, amap)
        if len(amap_ref) >= MIN_COMMON_ATOMS:
            p2, r2 = zip(*amap_ref)
            rms2, _ = _kabsch(pred_xyz[list(p2)], ref_xyz[list(r2)])
            if rms2 < rmsd_full:
                rmsd_refined = rms2
                n_common_refined = len(amap_ref)

    return {
        "rmsd_full": rmsd_full,
        "n_common": len(amap),
        "n_ref": len(ref_xyz),
        "rmsd_refined": rmsd_refined,
        "n_common_refined": n_common_refined,
    }


# ─────────────────────── PDBQT → PDB conversion ───────────────────────

def convert_pdbqt_to_clean_pdb(src: Path, dst: Path, template: Path | None = None) -> None:
    """
    Convert a PDBQT file to clean PDB, with optional bond-order template.
    """
    lines: list[str] = []
    for ln in src.read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            elem = _guess_element(ln[12:16])
            fixed = f"{ln[:54]}  1.00  0.00          {elem:>2}\n"
            lines.append(fixed)
    if not lines:
        raise ValueError("No ATOM/HETATM records – invalid PDBQT?")
    dst.write_text("".join(lines) + "END\n")

    # Apply template bond orders if provided
    if template and template.is_file():
        tmpl = Chem.MolFromPDBFile(str(template), removeHs=False)
        mol  = Chem.MolFromPDBFile(str(dst),      removeHs=False, proximityBonding=True)
        if tmpl and mol:
            try:
                mol = AllChem.AssignBondOrdersFromTemplate(tmpl, mol)
                Chem.SanitizeMol(mol, catchErrors=True)
                Chem.MolToPDBFile(mol, str(dst), flavor=0)
            except Exception:
                pass


# ─────────────────────── Overlay drawing ───────────────────────────────

def _draw_overlay_svg(ref: Chem.Mol, pred: Chem.Mol, fp: Path,
                      w: int = 450, h: int = 450, scale: float = 1.0) -> None:
    """
    Draw a 2D overlay of reference and predicted molecules to SVG.
    """
    drawer = rdMolDraw2D.MolDraw2DSVG(int(w*scale), int(h*scale))
    opts = drawer.drawOptions()
    opts.padding = 0.02
    opts.fixedBondLength = 22
    opts.bondLineWidth = 1.0
    opts.useBWAtomPalette()

    grey = {i: (0.3,0.3,0.3,0.25) for i in range(ref.GetNumAtoms())}
    red  = {i: (0.9,0.1,0.1,0.6) for i in range(pred.GetNumAtoms())}

    rdMolDraw2D.PrepareAndDrawMolecule(drawer, ref,
                                       highlightAtoms=list(grey),
                                       highlightAtomColors=grey,
                                       highlightBondColors=grey)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, pred,
                                       highlightAtoms=list(red),
                                       highlightAtomColors=red,
                                       highlightBondColors=red)
    drawer.FinishDrawing()
    fp.write_text(drawer.GetDrawingText())


def _find_native_file(native_dir: Path, rec_id: str, lig_id: str) -> Path | None:
    """
    Locate native ligand file: prefer .pdb, fallback to .pdbqt.
    """
    pdb = native_dir / f"{rec_id}_{lig_id}.pdb"
    if pdb.exists():
        return pdb
    pdbqt = native_dir / f"{rec_id}_{lig_id}.pdbqt"
    return pdbqt if pdbqt.exists() else None


# ─────────────────────── Main API ───────────────────────────────────────

def run_rmsd_and_plot(cfg: dict, log: logging.Logger) -> None:
    """
    Main entry point: compute RMSD and generate overlays for all
    native_redock PDBQT files under output_folder.

    Creates per-complex directories under `visuals`, writes overlays
    and a summary CSV with RMSD metrics.
    """
    native_dir = Path(cfg["paths"]["native_ligands_folder"])
    out_dir    = Path(cfg["paths"]["output_folder"])
    vis_root   = Path(cfg["paths"]["visuals"])
    vis_root.mkdir(parents=True, exist_ok=True)

    # Find all native-redock poses
    dock_files = list(out_dir.rglob("*__*__native_redock.pdbqt"))
    if not dock_files:
        log.error("No *_native_redock.pdbqt files found under %s", out_dir)
        return
    log.info("Found %d native_redock poses", len(dock_files))

    pattern = re.compile(r"(.+?)__(.+?)__native_redock\.pdbqt$", re.IGNORECASE)
    rows: list[dict] = []

    for dock in dock_files:
        m = pattern.match(dock.name)
        if not m:
            log.debug("Skip (name mismatch): %s", dock.name)
            continue
        rec_id, lig_id = m.groups()
        stem    = f"{rec_id}__{lig_id}__native_redock"
        vis_dir = vis_root / f"native_redock__{rec_id}__{lig_id}"
        vis_dir.mkdir(parents=True, exist_ok=True)

        native_fp = _find_native_file(native_dir, rec_id, lig_id)
        if not native_fp:
            log.warning("Missing native ligand for %s_%s", rec_id, lig_id)
            continue

        # 1) Compute RMSD metrics
        try:
            metrics = calc_rmsd_from_pdbqt(native_fp, dock, heavy_only=True)
            msg = f"RMSD={metrics['rmsd_full']:.3f} Å ({metrics['n_common']}/{metrics['n_ref']})"
            if metrics['rmsd_refined'] is not None:
                msg += f" | refined={metrics['rmsd_refined']:.3f} Å ({metrics['n_common_refined']} atoms)"
            log.info("%s: %s", stem, msg)
        except Exception as exc:
            log.error("RMSD error for %s: %s", dock.name, exc)
            continue

        # 2) Generate 2D overlay
        try:
            # Convert native to PDB if needed
            if native_fp.suffix.lower() == ".pdb":
                ref_pdb = native_fp
            else:
                ref_pdb = vis_dir / f"{rec_id}_{lig_id}_native.pdb"
                convert_pdbqt_to_clean_pdb(native_fp, ref_pdb, template=None)

            pred_pdb = dock.with_suffix(".pdb")
            convert_pdbqt_to_clean_pdb(dock, pred_pdb,
                                       template=ref_pdb if ref_pdb.exists() else None)

            ref_m  = Chem.MolFromPDBFile(str(ref_pdb),  removeHs=False, sanitize=False)
            pred_m = Chem.MolFromPDBFile(str(pred_pdb), removeHs=False, sanitize=False)
            if ref_m and pred_m:
                svg_fp = vis_dir / f"{stem}_overlay.svg"
                png_fp = vis_dir / f"{stem}_overlay.png"
                pdf_fp = vis_dir / f"{stem}_overlay.pdf"
                _draw_overlay_svg(ref_m, pred_m, svg_fp)
                cairosvg.svg2png(url=str(svg_fp), write_to=str(png_fp))
                cairosvg.svg2pdf(url=str(svg_fp), write_to=str(pdf_fp))
                log.info("Overlay saved for %s", stem)
            else:
                log.warning("Overlay skipped for %s (RDKit parse fail)", stem)
        except Exception as exc:
            log.warning("Overlay skipped for %s: %s", stem, exc)

        rows.append({
            "Complex": stem,
            "RMSD": metrics['rmsd_full'],
            "n_common": metrics['n_common'],
            "n_ref": metrics['n_ref'],
            "RMSD_refined": metrics['rmsd_refined'],
            "n_common_refined": metrics['n_common_refined'],
        })

    # Write summary CSV
    if rows:
        df = pd.DataFrame(rows).sort_values("RMSD")
        summary_fp = vis_root / "rmsd_summary.csv"
        df.to_csv(summary_fp, index=False)
        high = df[df.RMSD > 2.0]
        if not high.empty:
            log.info("HIGH RMSD (>2Å): " +
                     ", ".join(f"{c}:{r:.2f}" for c,r in zip(high.Complex, high.RMSD)))
        log.info("RMSD summary saved: %s", summary_fp)
    else:
        log.info("No RMSD records generated.")




