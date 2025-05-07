from pathlib import Path
import subprocess, tempfile, shutil
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

# ── altLoc splitter ────────────────────────────────────
def strip_alt_conformations_with_adt(pdb: Path) -> list[Path]:
    """Run ADT script to split altLoc into *_A.pdb, *_B.pdb, etc. Returns list of generated files."""
    subprocess.run(['python', str(split_alt_conf), '-r', str(pdb)], check=True)
    generated = list(pdb.parent.glob(f"{pdb.stem}_?.pdb"))
    if not generated:
        raise FileNotFoundError(f"No altLoc-split files found for {pdb.name}")
    return generated

def _strip_hetatm(pdb_path: Path) -> Path:
    out = pdb_path.with_name(pdb_path.stem + "_clean.pdb")
    with open(pdb_path) as f_in, open(out, 'w') as f_out:
        for line in f_in:
            if line.startswith("ATOM"):
                f_out.write(line)
    return out

def _fix_pdbqt_format(pdbqt_path: Path):
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

# ── główna logika ──────────────────────────────────────
def select_and_clean_best_conformation(pdb_id: str, raw_dir: Path, cleaned_dir: Path) -> Path:
    """
    Wybiera najlepszą konformację (preferując *_A.pdb), czyści HETATM i zapisuje do cleaned_dir/{pdb_id}.pdb.
    """
    candidates = list(raw_dir.glob(f"{pdb_id}_?.pdb"))
    if not candidates:
        raise FileNotFoundError(f"No altLoc-split files found for {pdb_id} in {raw_dir}")

    chain_a = [p for p in candidates if "_A" in p.stem]
    selected = chain_a[0] if chain_a else candidates[0]

    cleaned_dir.mkdir(exist_ok=True, parents=True)
    stripped = _strip_hetatm(selected)
    final_path = cleaned_dir / f"{pdb_id}.pdb"
    shutil.move(str(stripped), final_path)

    return final_path

def receptor_to_pdbqt(pdb: Path, out_dir: Path, strip_waters: bool = True):
    """Convert receptor PDB to PDBQT using AutoDockTools. Assumes altLoc is already resolved."""
    out = out_dir / f"{pdb.stem}.pdbqt"
    cmd = [
        'python', str(prep_receptor),
        '-r', str(pdb),
        '-o', str(out),
        '-A', 'checkhydrogens',
        '-U', 'nphs_lps_waters_nonstdres'
    ]
    if strip_waters:
        cmd += ['-U', 'waters']

    subprocess.run(cmd, check=True)
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

