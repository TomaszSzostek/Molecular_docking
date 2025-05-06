"""Download PDB files by ID list (requires requests)."""
import requests, gzip, shutil
from pathlib import Path


def fetch(pdb_id: str, out_dir: Path):
    url = f'https://files.rcsb.org/download/{pdb_id.upper()}.pdb.gz'
    gz_path = out_dir / f'{pdb_id}.pdb.gz'
    pdb_path = out_dir / f'{pdb_id}.pdb'
    r = requests.get(url, timeout=30)
    gz_path.write_bytes(r.content)
    with gzip.open(str(gz_path),'rb') as gzf, open(str(pdb_path),'wb') as out:
        shutil.copyfileobj(gzf,out)
    gz_path.unlink()
    return pdb_path

from Bio.PDB import PDBParser, PDBIO, Select

class LigandSelect(Select):
    def __init__(self, het_id):
        self.het_id = het_id
    def accept_residue(self, residue):
        return residue.id[0].strip()=="H_" and residue.resname==self.het_id

class ProteinSelect(Select):
    def accept_residue(self, residue):
        return residue.id[0]==" "  # standard AA only


def split_receptor_ligand(pdb_path: Path, ligand_resname: str, out_receptor: Path, out_ligand: Path):
    struct = PDBParser(QUIET=True).get_structure('rec', str(pdb_path))
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(out_ligand), LigandSelect(ligand_resname))
    io.save(str(out_receptor), ProteinSelect())