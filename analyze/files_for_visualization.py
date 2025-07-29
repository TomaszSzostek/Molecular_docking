"""
files_for_visualization.py

Convert docking output PDBQTs into organized PDB complexes and PLIP reports
for downstream visualization in py3Dmol.

This module provides helper routines to:
  1. Extract ATOM/HETATM records and generate clean PDB files.
  2. Assign unique atom serials and residue metadata to ligand PDBs.
  3. Assemble receptor and ligand into a single complex file.
  4. Produce PLIP-compatible PDB inputs and invoke PLIP for interaction analysis.

All routines preserve existing file structure and skip already-generated artifacts.
"""

from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET
import logging
from openbabel import pybel

log = logging.getLogger(__name__)

# Color codes for PLIP interactions (used by external viewers)
COLOR = {
    "hb": "red",     # hydrogen bond
    "hp": "yellow",  # hydrophobic contact
    "pc": "purple",  # pi–cation
    "ps": "orange",  # pi–stacking
    "sb": "blue",    # salt bridge
    "wb": "cyan",    # water-mediated interaction
    "ha": "lime",    # halogen bond
    "me": "orange",  # metal coordination (edge)
    "mc": "orange",  # metal coordination (center)
}

# ────────────────────────────────────────────────────────── helpers ──────────────────────────────────────────────────────────

def _pdbqt_to_pdb(src: Path, dst: Path):
    """
    Copy ATOM/HETATM lines from a PDBQT file and append TER/END records.

    Parameters
    ----------
    src : Path
        Source .pdbqt file with AutoDock atom-type annotations.
    dst : Path
        Destination .pdb file that will be created or overwritten.
    """
    with src.open() as fin, dst.open("w") as fout:
        for ln in fin:
            if ln.startswith(("ATOM", "HETATM")):
                # retain columns 1-54 (coordinates & metadata) and strip ADT types
                fout.write(f"{ln[:54].rstrip():<54}\n")
        # close the PDB with TER and END markers
        fout.write("TER\nEND\n")


def _max_serial(pdb_path: Path) -> int:
    """
    Return the highest atom serial number in a PDB file.

    Parameters
    ----------
    pdb_path : Path
        Path to a .pdb file to scan for atom serials (cols 7-11).

    Returns
    -------
    int
        Maximum serial found, or 0 if none.
    """
    last = 0
    with pdb_path.open() as fh:
        for ln in fh:
            if ln.startswith(("ATOM", "HETATM")):
                try:
                    serial = int(ln[6:11])
                    last = max(last, serial)
                except ValueError:
                    # skip lines with non-numeric serial fields
                    pass
    return last


def _write_ligand_pdb(
    lig_pdbqt: Path,
    dst: Path,
    start_serial: int,
    resname: str = "LIG",
    chain: str | None = None,
    resnum: int = 501,
    receptor_pdb: Path | None = None,
) -> None:
    """
    Generate a PDB from a PDBQT ligand, renumber atoms, and set residue metadata.

    Parameters
    ----------
    lig_pdbqt : Path
        Ligand file in .pdbqt format.
    dst : Path
        Output .pdb file path.
    start_serial : int
        Starting atom serial number to ensure uniqueness.
    resname : str, optional
        Three-letter residue name (default 'LIG').
    chain : str or None, optional
        Chain identifier; if None or in use, next free A-Z is chosen.
    resnum : int, optional
        Residue sequence number (default 501).
    receptor_pdb : Path or None, optional
        PDB of receptor to detect used chains; if None, assumes '<dst.parent>/receptor.pdb'.
    """
    import tempfile, os

    # determine which chains are occupied in the receptor
    if receptor_pdb is None:
        receptor_pdb = dst.parent / "receptor.pdb"

    used_chains: set[str] = set()
    if receptor_pdb.is_file():
        with receptor_pdb.open() as fh:
            used_chains = {ln[21] for ln in fh if ln.startswith(("ATOM", "HETATM"))}

    # pick a free chain letter if needed
    if chain in (None, "", "None") or chain in used_chains:
        chain = next((c for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" if c not in used_chains), "Z")

    # use OpenBabel to read PDBQT and write a PDB with residue info
    tmp = Path(tempfile.mkstemp(suffix=".pdb")[1])
    mol = next(pybel.readfile("pdbqt", str(lig_pdbqt)))

    for atom in mol:
        res = atom.OBAtom.GetResidue() or mol.OBMol.NewResidue()
        res.SetName(resname)
        res.SetChain(chain)
        res.SetNum(resnum)
        atom.OBAtom.SetResidue(res)
    mol.OBMol.SetTitle(resname)
    mol.write("pdb", str(tmp), overwrite=True)

    # renumber atom serials and write final PDB
    serial = start_serial
    with tmp.open() as fin, dst.open("w") as fout:
        for ln in fin:
            if ln.startswith(("ATOM", "HETATM")):
                ln = f"{ln[:6]}{serial:>5}{ln[11:]}"
                serial += 1
            fout.write(ln)

    os.remove(tmp)  # cleanup temporary file


# ───────────── prepare_complex_pdb ─────────────

def prepare_complex_pdb(rec_pdbqt: Path, lig_pdbqt: Path, complex_dir: Path):
    """
    Assemble receptor and ligand PDBs into a visualization complex and PLIP input.

    Parameters
    ----------
    rec_pdbqt : Path
        Receptor file in .pdbqt format.
    lig_pdbqt : Path
        Ligand pose file '<rec>__<lig>__<tag>.pdbqt'.
    complex_dir : Path
        Directory where output files will be created.

    Returns
    -------
    tuple[Path, Path, Path, Path]
        Paths to receptor.pdb, ligand.pdb, plip_pdb, view_pdb files.
    """
    rec_id = rec_pdbqt.stem
    lig_id, tag = lig_pdbqt.stem.split("__")[1:3]

    rec_pdb  = complex_dir / f"{rec_id}_receptor.pdb"
    lig_pdb  = complex_dir / f"{rec_id}_{lig_id}_ligand.pdb"
    plip_pdb = complex_dir / f"{rec_id}_{lig_id}_plip.pdb"
    view_pdb = complex_dir / f"complex__{rec_id}__{lig_id}__{tag}.pdb"

    # 1) create receptor PDB if missing
    if not rec_pdb.exists():
        _pdbqt_to_pdb(rec_pdbqt, rec_pdb)

    # 2) create ligand PDB with renumbering
    if not lig_pdb.exists():
        offset = _max_serial(rec_pdb) + 1
        _write_ligand_pdb(lig_pdbqt, lig_pdb, start_serial=offset)

    # 3) prepare PLIP input by concatenating receptor + ligand
    if not plip_pdb.exists():
        with plip_pdb.open("w") as fout:
            # write receptor atoms
            for ln in rec_pdb.read_text().splitlines():
                if ln.startswith(("ATOM", "HETATM")):
                    fout.write(ln + "\n")
            fout.write("TER\n")
            # write ligand atoms
            for ln in lig_pdb.read_text().splitlines():
                if ln.startswith(("ATOM", "HETATM")):
                    fout.write(ln + "\n")
            fout.write("END\n")

    # 4) prepare simplified PDB for py3Dmol visualization
    if not view_pdb.exists():
        with rec_pdb.open() as frec, lig_pdb.open() as flig, view_pdb.open("w") as fout:
            for ln in frec:
                if ln.startswith(("ATOM", "HETATM")):
                    fout.write(ln.rstrip() + "\n")
            fout.write("TER\n")
            for ln in flig:
                if ln.startswith(("ATOM", "HETATM")):
                    fout.write(ln.rstrip() + "\n")
            fout.write("END\n")

    return rec_pdb, lig_pdb, plip_pdb, view_pdb


def generate_files(cfg: dict, log):
    """
    Discover docking results and generate visualization artifacts.

    Steps:
      1) Identify PDBQT ligand results matching the docking mode.
      2) For each, ensure output directory exists.
      3) Call prepare_complex_pdb to build PDBs.
      4) Invoke PLIP to create interaction XML reports.

    Parameters
    ----------
    cfg : dict
        Pipeline configuration with 'paths' and 'docking_mode'.
    log : Logger
        Logger for progress and errors.

    Returns
    -------
    None
    """
    viz_cfg = cfg.get("visualization", {})
    if not viz_cfg.get("enabled", True):
        log.info("Visualization disabled – skipping")
        return

    mode = cfg.get("docking_mode", "matrix")
    out_dir = Path(cfg["paths"]["output_folder"])
    base_vis = Path(cfg["paths"]["visuals"])
    rec_dir = Path(cfg["paths"]["receptors_cleaned_folder"])

    if mode == "redock_native":
        result_dir = out_dir / "dock_native"
        expected_tag = "native_redock"
    elif mode == "diagonal":
        result_dir = out_dir / "diagonal"
        expected_tag = "diag"
    else:
        result_dir = out_dir / "matrix"
        expected_tag = "dock"

    lig_files = list(result_dir.rglob("*__*__*.pdbqt"))
    if not lig_files:
        log.warning("No docking results (*.pdbqt) for mode %r – skipping", mode)
        return

    log.info("Preparing complexes for %d ligands (mode=%s)…", len(lig_files), mode)

    for lig_pdbqt in lig_files:
        rec_id, lig_id, tag = lig_pdbqt.stem.split("__")[:3]
        if tag != expected_tag:
            continue

        rec_pdbqt = rec_dir / f"{rec_id}.pdbqt"
        if not rec_pdbqt.exists():
            log.warning("Receptor %s missing – skipping %s", rec_id, lig_pdbqt.name)
            continue

        complex_dir = base_vis / f"{tag}__{rec_id}__{lig_id}"
        complex_dir.mkdir(parents=True, exist_ok=True)

        out_pdb = complex_dir / f"complex__{rec_id}__{lig_id}__{tag}.pdb"
        if not out_pdb.exists():
            try:
                prepare_complex_pdb(rec_pdbqt, lig_pdbqt, complex_dir)
                log.debug("Prepared complex: %s", out_pdb.name)
            except Exception as e:
                log.error("Failed to prepare complex for %s: %s", lig_pdbqt.name, e)
                continue

        plip_dir = complex_dir / "plip"
        plip_xml = plip_dir / "report.xml"
        if not plip_xml.exists():
            plip_dir.mkdir(exist_ok=True)
            try:
                subprocess.run([
                    "plip", "-f", str(out_pdb), "-o", str(plip_dir), "-x"
                ], check=True)
                log.debug("Generated PLIP XML for: %s", out_pdb.name)
            except Exception as e:
                log.error("PLIP failed for %s: %s", out_pdb.name, e)








