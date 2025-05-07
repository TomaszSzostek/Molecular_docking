"""Convert a CSV of SMILES → PDBQT and save to ligands folder."""
from pathlib import Path
import pandas as pd, os, tempfile
from rdkit import Chem
from .adt_preparator import ligand_to_pdbqt

def csv_to_pdbqt(csv_path: Path, smiles_col: str, id_col: str, out_dir: Path):
    df = pd.read_csv(csv_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row[smiles_col]))
        if mol is None:
            continue
        Chem.AddHs(mol)
        tmp = Path(tempfile.mkstemp(suffix=".pdb")[1])
        Chem.MolToPDBFile(mol, str(tmp))
        ligand_to_pdbqt(tmp, out_dir, sanitize=True)
        os.remove(tmp)
