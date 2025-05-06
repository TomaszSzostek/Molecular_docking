from pathlib import Path
import subprocess, shutil

# locate ADT scripts inside env
from importlib.util import find_spec
adt_path = Path(find_spec('AutoDockTools_py3').origin).parent
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