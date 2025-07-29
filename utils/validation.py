from pathlib import Path
import numpy as np

def files_ok(lig_path: Path) -> bool:
    """Check if docking output files exist for a ligand."""
    stem = lig_path.stem
    log_file = lig_path.with_suffix(".log")
    pdbqt_out = lig_path  # already ends with .pdbqt
    return pdbqt_out.exists() and log_file.exists() and pdbqt_out.stat().st_size > 0


def rmsd_atoms(ref_pdb: str, pose_pdb: str) -> float:
    """Very simple heavy‑atom RMSD (BioPython-free to keep deps light)."""
    def coords(pdb):
        return np.array([[float(p[x:x+8]) for x in (30,38,46)]
                         for p in open(pdb) if p.startswith("ATOM")])
    a, b = coords(ref_pdb), coords(pose_pdb)
    if len(a) != len(b):
        return np.nan
    return np.sqrt(((a-b)**2).sum(axis=1).mean())
