"""Generate grid box JSON (center,size) from native ligand COM."""
import json, numpy as np
from pathlib import Path
from rdkit import Chem


def com_xyz(pdbqt: Path):
    mol = Chem.MolFromPDBFile(str(pdbqt), removeHs=False)
    conf = mol.GetConformer()
    coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    return coords.mean(axis=0)


def make_box(receptor_pdbqt: Path, native_lig_pdbqt: Path, margin: float=4.0):
    center = com_xyz(native_lig_pdbqt)
    size = np.array([margin,margin,margin])*2  # cubic box
    meta = {'center': center.tolist(), 'size': size.tolist()}
    with open(receptor_pdbqt.with_suffix('.box.json'),'w') as fh:
        json.dump(meta,fh,indent=2)
    return meta