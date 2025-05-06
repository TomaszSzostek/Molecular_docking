from pathlib import Path
import subprocess, shutil, os
from pathlib import Path
import subprocess, tempfile
from importlib.util import find_spec
from rdkit import Chem
from rdkit.Chem import AllChem

# locate ADT scripts inside env
adt_path = Path(find_spec('AutoDockTools_py3').origin).parent
prep_receptor = adt_path / 'prepare_receptor4.py'
prep_ligand   = adt_path / 'prepare_ligand4.py'

# ── Receptor ─────────────────────────────────────────

def receptor_to_pdbqt(pdb: Path, out_dir: Path, strip_waters: bool = True):
    """Convert PDB → PDBQT using ADT.  If *strip_waters* remove HOH & add hydrogens."""
    out = out_dir / f"{pdb.stem}.pdbqt"
    cmd = ['python', str(prep_receptor), '-r', str(pdb), '-o', str(out), '-A', 'hydrogens']
    if strip_waters:
        cmd += ['-U', 'waters']
    subprocess.run(cmd, check=True)
    return out

# ── Ligand ───────────────────────────────────────────

def _sanitize_ligand(pdb_in: Path) -> Path:
    """RDKit sanitize & UFF minimise; returns path to temp PDB."""
    mol = Chem.MolFromPDBFile(str(pdb_in), removeHs=False)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    tmp = Path(tempfile.mkstemp(suffix='.pdb')[1])
    Chem.MolToPDBFile(mol, str(tmp))
    return tmp


def ligand_to_pdbqt(pdb: Path, out_dir: Path, sanitize: bool = True):
    """Convert PDB → PDBQT via ADT; optional RDKit pre‑sanitise."""
    src_pdb = _sanitize_ligand(pdb) if sanitize else pdb
    out = out_dir / f"{pdb.stem}.pdbqt"
    subprocess.run(['python', str(prep_ligand), '-l', str(src_pdb), '-o', str(out), '-A', 'hydrogens'], check=True)
    if sanitize:
        src_pdb.unlink(missing_ok=True)
    return out

"""Wrap AutoDockTools `prepare_receptor4.py` / `prepare_ligand4.py`."""
from pathlib import Path
import subprocess, shutil

# locate ADT scripts inside env
from importlib.util import find_spec
_custom = os.getenv('ADT_HOME')
if _custom and Path(_custom).is_dir():
    adt_path = Path(_custom)
else:
    try:
        from importlib.util import find_spec
        adt_path = Path(find_spec('AutoDockTools_py3').origin).parent
    except (ModuleNotFoundError, AttributeError):
        raise RuntimeError("AutoDockTools scripts not found. Install the package or set ADT_HOME env var to its directory.")

prep_receptor = adt_path / 'prepare_receptor4.py'
prep_ligand   = adt_path / 'prepare_ligand4.py'


def receptor_to_pdbqt(pdb: Path, out_dir: Path):
    out = out_dir / f"{pdb.stem}.pdbqt"
    subprocess.run(['python', str(prep_receptor), '-r', str(pdb), '-o', str(out), '-A', 'hydrogens'], check=True)
    return out


def ligand_to_pdbqt(pdb: Path, out_dir: Path):
    out = out_dir / f"{pdb.stem}.pdbqt"
    subprocess.run(['python', str(prep_ligand), '-l', str(pdb), '-o', str(out), '-A', 'hydrogens'], check=True)
    return out