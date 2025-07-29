"""
prepare_inputs/adt_preparator.py

Wrap AutoDockTools Py3 scripts to prepare receptor and ligand PDBQT files.

This module provides:
  - strip_alt_conformations_with_adt: remove alternate locations via ADT script.
  - _strip_hetatm: filter HETATM records, waters and ions based on proximity.
  - _fix_pdbqt_format: clean malformed ATOM/HETATM lines in PDBQT.
  - receptor_to_pdbqt: full receptor preparation pipeline calling ADT.
  - _sanitize and ligand_to_pdbqt: optional RDKit sanitization and ligand PDBQT conversion.

Public API:
    receptor_to_pdbqt(pdb: Path, out_dir: Path, strip_waters: bool=True,
                      strip_hetatm: bool=True, ligand_pdb: Optional[Path]=None) -> Path
    ligand_to_pdbqt(pdb: Path, out_dir: Path, sanitize: bool=False) -> Path
"""

from pathlib import Path
import subprocess, tempfile
from rdkit import Chem
from rdkit.Chem import AllChem
from importlib.util import find_spec
from typing import Optional

# ── locate ADT prepare_inputs scripts from installed package ──
spec = find_spec("AutoDockTools_py3") or find_spec("AutoDockTools")
if spec is None:
    raise RuntimeError("AutoDockTools not found in current Python environment.")

adt_root = Path(spec.origin).parent
prep_receptor = adt_root / "Utilities24" / "prepare_receptor4.py"
prep_ligand   = adt_root / "Utilities24" / "prepare_ligand4.py"
split_alt_conf = adt_root / "Utilities24" / "prepare_pdb_split_alt_confs.py"


def strip_alt_conformations_with_adt(pdb: Path, prefer=("A", "")) -> Path:
    """
    Run ADT split-alt-locs and keep only the preferred alternate conformation.

    Parameters
    ----------
    pdb : Path
        Input .pdb file with potential alternate locations.
    prefer : tuple of str
        Preferred altLoc tags in order (default ('A', '')).

    Returns
    -------
    Path
        The same Path pointing to the file with alternates stripped.
    """
    # 1. Run ADT to generate files like mypdb_A.pdb, mypdb_B.pdb
    subprocess.run(['python', str(split_alt_conf), '-r', str(pdb)], check=True)

    # 2. Gather all alternate-location files
    alt_files = list(pdb.parent.glob(f"{pdb.stem}_?.pdb"))
    if not alt_files:
        # no alternates → return original
        return pdb

    # 3. Choose preferred altLoc file
    chosen: Path | None = None
    for tag in prefer:
        candidate = pdb.with_name(f"{pdb.stem}_{tag or ' '}.pdb")
        if candidate.exists():
            chosen = candidate
            break
    if chosen is None:
        chosen = alt_files[0]

    # 4. Replace original with chosen and remove others
    target = pdb
    chosen.replace(target)
    for f in alt_files:
        if f != chosen:
            f.unlink(missing_ok=True)

    return target


# ──────────────────────────────────────────────────────────────────────
# 1.  _strip_hetatm
# ──────────────────────────────────────────────────────────────────────

def _strip_hetatm(
    pdb_path: Path,
    ligand_pdb: Optional[Path] = None,
    cutoff: float = 3.0
) -> Path:
    """
    Remove unwanted HETATM records, keep waters/ions near ligand if provided.

    Parameters
    ----------
    pdb_path : Path
        Input receptor PDB file.
    ligand_pdb : Optional[Path]
        Ligand PDB for proximity filtering.
    cutoff : float
        Distance cutoff in Å for keeping water/ions (default 3.0).

    Returns
    -------
    Path
        Path to cleaned PDB with suffix '_clean.pdb'.
    """
    # 1) Collect ligand atom coordinates if ligand provided
    lig_coords: list[tuple[float, float, float]] = []
    if ligand_pdb is not None and ligand_pdb.exists():
        with ligand_pdb.open() as lf:
            for ln in lf:
                if ln.startswith(("ATOM", "HETATM")) and len(ln) >= 54:
                    try:
                        lig_coords.append((
                            float(ln[30:38]),
                            float(ln[38:46]),
                            float(ln[46:54])
                        ))
                    except ValueError:
                        pass
    cutoff2 = cutoff * cutoff  # squared cutoff

    # 2) Open input and output files
    out = pdb_path.with_name(pdb_path.stem + "_clean.pdb")
    water_tags = {"HOH", "WAT", "DOD"}
    special    = {"ZN", "MG"}  # metal ions to consider by distance

    with pdb_path.open() as fin, out.open("w") as fout:
        for ln in fin:
            if len(ln) < 54:
                # skip non-ATOM/HETATM lines
                continue

            rec = ln[:6]
            alt = ln[16]
            res = ln[17:20].strip().upper()
            keep_atom = False

            # --- ATOM records ---
            if rec == "ATOM  " and alt in (" ", "A"):
                keep_atom = True

            # --- HETATM records ---
            elif rec == "HETATM" and alt in (" ", "A"):

                # waters and special ions
                if res in water_tags | special:
                    if lig_coords:
                        # keep if within cutoff of any ligand atom
                        try:
                            x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                            keep_atom = any(
                                (x - lx)**2 + (y - ly)**2 + (z - lz)**2 <= cutoff2
                                for lx, ly, lz in lig_coords
                            )
                        except ValueError:
                            pass
                    elif res in special:
                        # no ligand provided: keep Zn/Mg only
                        keep_atom = True

                # all other HETATM except waters
                elif res not in water_tags:
                    keep_atom = True

            if keep_atom:
                fout.write(ln)

    return out


def _fix_pdbqt_format(pdbqt_path: Path) -> None:
    """
    Remove malformed ATOM/HETATM lines from a PDBQT file.

    Parameters
    ----------
    pdbqt_path : Path
        Path to the .pdbqt file to clean.
    """
    valid_lines = []
    with open(pdbqt_path, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM", "ROOT", "BRANCH", "ENDROOT", "TORSDOF", "REMARK")):
                try:
                    # validate floating-point columns
                    float(line[30:38]); float(line[38:46]); float(line[46:54])
                    valid_lines.append(line)
                except ValueError:
                    # skip invalid coordinate lines
                    continue
    with open(pdbqt_path, 'w') as f:
        f.writelines(valid_lines)


# ──────────────────────────────────────────────────────────────────────
# 2.  receptor_to_pdbqt
# ──────────────────────────────────────────────────────────────────────

def receptor_to_pdbqt(
    pdb: Path,
    out_dir: Path,
    strip_waters: bool = True,
    strip_hetatm: bool = True,
    ligand_pdb: Optional[Path] = None
) -> Path:
    """
    Convert a receptor PDB file to PDBQT using AutoDockTools.

    Steps:
      1. Strip alternate conformations.
      2. Clean HETATM/waters near ligand if requested.
      3. Run ADT prepare_receptor4.py script.
      4. Remove temporary clean PDB.
      5. Fix malformed lines in final PDBQT.

    Parameters
    ----------
    pdb : Path
        Input receptor PDB file.
    out_dir : Path
        Directory to write .pdbqt file.
    strip_waters : bool
        If True, remove water molecules via ADT '-U waters'.
    strip_hetatm : bool
        If True, apply _strip_hetatm filtering.
    ligand_pdb : Optional[Path]
        Ligand PDB for proximity-based HETATM filtering.

    Returns
    -------
    Path
        Path to generated .pdbqt file.
    """
    # 1) select single altLoc
    pdb_cleaned = strip_alt_conformations_with_adt(pdb)

    # 2) filter waters/HETATM
    pdb_in = (_strip_hetatm(pdb_cleaned, ligand_pdb=ligand_pdb)
              if strip_hetatm else pdb_cleaned)

    out = out_dir / f"{pdb.stem}.pdbqt"

    # 3) call ADT prepare_receptor4.py
    cmd = [
        'python', str(prep_receptor),
        '-r', str(pdb_in),
        '-o', str(out),
        '-A', 'checkhydrogens'
    ]
    if strip_waters:
        cmd += ['-U', 'waters']
    cmd += ['-U', 'nphs_lps']

    subprocess.run(cmd, check=True)

    # 4) remove temporary clean PDB
    if strip_hetatm:
        pdb_in.unlink(missing_ok=True)

    # 5) correct malformed lines
    _fix_pdbqt_format(out)
    return out


# ── ligand conversion ───────────────────────────────

def _sanitize(pdb_in: Path) -> Path:
    """
    Sanitize and embed a ligand PDB using RDKit, returns a temporary PDB.

    Parameters
    ----------
    pdb_in : Path
        Input ligand PDB file.

    Returns
    -------
    Path
        Temporary sanitized PDB path.
    """
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


def ligand_to_pdbqt(pdb: Path, out_dir: Path, sanitize: bool = False) -> Path:
    """
    Convert a ligand PDB to PDBQT using AutoDockTools, with optional RDKit sanitization.

    Parameters
    ----------
    pdb : Path
        Input ligand PDB file.
    out_dir : Path
        Directory to write the .pdbqt file.
    sanitize : bool
        If True, run RDKit sanitization and embedding before conversion.

    Returns
    -------
    Path
        Path to generated .pdbqt file.
    """
    src = _sanitize(pdb) if sanitize else pdb
    out = out_dir / f"{pdb.stem.replace('_receptor','')}.pdbqt"

    subprocess.run(["python", str(prep_ligand), "-l", str(src), "-o",
                    str(out), "-A", "hydrogens",
                    "no_tors"], check=True)
    if sanitize:
        src.unlink(missing_ok=True)

    return out
