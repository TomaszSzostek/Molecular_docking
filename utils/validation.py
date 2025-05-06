from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

REQUIRED = ['poses.pdbqt','log.txt','atom_terms.txt']

def files_ok(directory: Path):
    for fname in REQUIRED:
        f = directory / fname
        if not f.is_file() or f.stat().st_size==0:
            return False
    return True


def rmsd(native_sdf: Path, docked_sdf: Path)->float:
    m1 = Chem.SDMolSupplier(str(native_sdf))[0]
    m2 = Chem.SDMolSupplier(str(docked_sdf))[0]
    return rdMolAlign.GetBestRMS(m1,m2)