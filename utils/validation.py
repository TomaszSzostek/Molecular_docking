from pathlib import Path
import numpy as np

def files_ok(dir_: Path):
    req = ["poses.pdbqt", "log.txt", "atom_terms.txt"]
    return all((dir_/f).is_file() and (dir_/f).stat().st_size > 0 for f in req)

def rmsd_atoms(ref_pdb: str, pose_pdb: str) -> float:
    """Very simple heavy‑atom RMSD (BioPython-free to keep deps light)."""
    def coords(pdb):
        return np.array([[float(p[x:x+8]) for x in (30,38,46)]
                         for p in open(pdb) if p.startswith("ATOM")])
    a, b = coords(ref_pdb), coords(pose_pdb)
    if len(a) != len(b):
        return np.nan
    return np.sqrt(((a-b)**2).sum(axis=1).mean())
