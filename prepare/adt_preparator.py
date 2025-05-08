from pathlib import Path
import subprocess, tempfile
from rdkit import Chem
from rdkit.Chem import AllChem
from importlib.util import find_spec

# ── locate ADT prepare scripts from installed package ──
spec = find_spec("AutoDockTools_py3") or find_spec("AutoDockTools")
if spec is None:
    raise RuntimeError("AutoDockTools not found in current Python environment.")

adt_root = Path(spec.origin).parent
prep_receptor = adt_root / "Utilities24" / "prepare_receptor4.py"
prep_ligand   = adt_root / "Utilities24" / "prepare_ligand4.py"
split_alt_conf = adt_root / "Utilities24" / "prepare_pdb_split_alt_confs.py"


def strip_alt_conformations_with_adt(pdb: Path) -> Path:
    """Remove alternate conformations using prepare_pdb_split_alt_confs.py, keeping only 'A'."""
    # Run ADT script to generate multiple files (one per altLoc)
    subprocess.run(['python', str(split_alt_conf), '-r', str(pdb)], check=True)

    # ADT script outputs file like: inputname_A.pdb
    alt_path = pdb.with_name(f"{pdb.stem}_A.pdb")

    if not alt_path.exists():
        raise FileNotFoundError(f"Expected output {alt_path.name} not found after altLoc splitting.")

    return alt_path


def _strip_hetatm(pdb_path: Path) -> Path:
    """Remove all non-ATOM lines and keep only altLoc=' ' or 'A'."""
    out = pdb_path.with_name(pdb_path.stem + "_clean.pdb")
    with open(pdb_path) as f_in, open(out, 'w') as f_out:
        for line in f_in:
            if line.startswith("ATOM") and (line[16] in (' ', 'A')):  # altLoc column
                f_out.write(line)
    return out


def _fix_pdbqt_format(pdbqt_path: Path):
    """Remove malformed ATOM lines with invalid coordinates."""
    valid_lines = []
    with open(pdbqt_path, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM", "ROOT", "BRANCH", "ENDROOT", "TORSDOF", "REMARK")):
                try:
                    float(line[30:38])
                    float(line[38:46])
                    float(line[46:54])
                    valid_lines.append(line)
                except ValueError:
                    continue
    with open(pdbqt_path, 'w') as f:
        f.writelines(valid_lines)


def receptor_to_pdbqt(pdb: Path, out_dir: Path, strip_waters: bool = True, strip_hetatm: bool = True):
    """Convert receptor PDB file to PDBQT format using AutoDockTools."""

    # 🧹 Usuwamy altLoc tylko z istniejącego finalnego pliku `1M17.pdb`
    pdb_cleaned = strip_alt_conformations_with_adt(pdb)
    pdb_in = _strip_hetatm(pdb_cleaned) if strip_hetatm else pdb_cleaned
    out = out_dir / f"{pdb.stem}.pdbqt"

    cmd = [
        'python', str(prep_receptor),
        '-r', str(pdb_in),
        '-o', str(out),
        '-A', 'checkhydrogens',
        '-U', 'nphs_lps_waters_nonstdres'
    ]
    if strip_waters:
        cmd += ['-U', 'waters']

    subprocess.run(cmd, check=True)

    if strip_hetatm:
        pdb_in.unlink(missing_ok=True)

    _fix_pdbqt_format(out)
    return out

# ── ligand conversion ───────────────────────────────
def _sanitize(pdb_in: Path) -> Path:
    mol = Chem.MolFromPDBFile(str(pdb_in), removeHs=False)
    if mol is None:
        raise ValueError(f"RDKit failed to parse ligand PDB: {pdb_in.name}")

    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    tmp = Path(tempfile.mkstemp(suffix=".pdb")[1])
    Chem.MolToPDBFile(mol, str(tmp))
    return tmp


def ligand_to_pdbqt(pdb: Path, out_dir: Path, sanitize: bool = True):
    src = _sanitize(pdb) if sanitize else pdb
    out = out_dir / f"{pdb.stem.replace('_receptor','')}.pdbqt"

    subprocess.run(["python", str(prep_ligand), "-l", str(src), "-o", str(out), "-A", "hydrogens"], check=True)
    if sanitize:
        src.unlink(missing_ok=True)

    return out
