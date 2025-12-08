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
import shutil
import xml.etree.ElementTree as ET
import logging
from typing import Optional
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


def _extract_waters_from_pdb(pdb_path: Path, ligand_coords: list[tuple[float, float, float]], 
                              cutoff: float = 5.0) -> list[str]:
    """
    Extract water molecules (HOH/WAT) from a PDB file that are within cutoff distance of ligand.
    
    Parameters
    ----------
    pdb_path : Path
        Path to original PDB file (may contain waters).
    ligand_coords : list[tuple[float, float, float]]
        List of (x, y, z) coordinates of ligand atoms.
    cutoff : float
        Maximum distance in Angstroms for water to be included (default 5.0).
    
    Returns
    -------
    list[str]
        List of PDB lines for water molecules within cutoff.
    """
    if not pdb_path.exists():
        return []
    
    waters = []
    cutoff2 = cutoff * cutoff
    water_resnames = {"HOH", "WAT", "DOD"}
    
    with pdb_path.open() as fh:
        for ln in fh:
            if not ln.startswith("HETATM"):
                continue
            resname = ln[17:20].strip().upper()
            if resname not in water_resnames:
                continue
            
            try:
                x = float(ln[30:38])
                y = float(ln[38:46])
                z = float(ln[46:54])
            except (ValueError, IndexError):
                continue
            
            # Check if water is within cutoff of any ligand atom
            for lx, ly, lz in ligand_coords:
                dist2 = (x - lx)**2 + (y - ly)**2 + (z - lz)**2
                if dist2 <= cutoff2:
                    waters.append(ln.rstrip())
                    break
    
    return waters


def _get_ligand_coords(lig_pdb: Path) -> list[tuple[float, float, float]]:
    """
    Extract coordinates of all ligand atoms from a PDB file.
    
    Parameters
    ----------
    lig_pdb : Path
        Path to ligand PDB file.
    
    Returns
    -------
    list[tuple[float, float, float]]
        List of (x, y, z) coordinates.
    """
    coords = []
    if not lig_pdb.exists():
        return coords
    
    with lig_pdb.open() as fh:
        for ln in fh:
            if ln.startswith(("ATOM", "HETATM")):
                try:
                    x = float(ln[30:38])
                    y = float(ln[38:46])
                    z = float(ln[46:54])
                    coords.append((x, y, z))
                except (ValueError, IndexError):
                    continue
    
    return coords


def _write_ligand_pdb(
    lig_pdbqt: Path,
    dst: Path,
    start_serial: int,
    resname: str = "LIG",
    chain: Optional[str] = None,
    resnum: int = 501,
    receptor_pdb: Optional[Path] = None,
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

def prepare_complex_pdb(rec_pdbqt: Path, lig_pdbqt: Path, complex_dir: Path, 
                        original_pdb: Optional[Path] = None):
    """
    Assemble receptor and ligand PDBs into a visualization complex and PLIP input.
    
    If original_pdb is provided, water molecules within 5Å of the ligand will be
    added to the PLIP input structure to enable water bridge detection.

    Parameters
    ----------
    rec_pdbqt : Path
        Receptor file in .pdbqt format.
    lig_pdbqt : Path
        Ligand pose file '<rec>__<lig>__<tag>.pdbqt'.
    complex_dir : Path
        Directory where output files will be created.
    original_pdb : Path or None, optional
        Original PDB file (may contain waters). If provided, waters near ligand
        will be extracted and added to PLIP input.

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

    # 3) prepare PLIP input by concatenating receptor + ligand + waters (if available)
    if not plip_pdb.exists():
        # Get ligand coordinates for water proximity filtering
        lig_coords = _get_ligand_coords(lig_pdb)
        
        # Extract waters from original PDB if available
        waters = []
        if original_pdb is not None and original_pdb.exists():
            waters = _extract_waters_from_pdb(original_pdb, lig_coords, cutoff=5.0)
            if waters:
                log.debug("Found %d water molecules near ligand in %s", len(waters), original_pdb.name)
        
        # Write PLIP input: receptor + ligand + waters
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
            # write waters if available
            if waters:
                fout.write("TER\n")
                for water_line in waters:
                    fout.write(water_line + "\n")
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
      1) Load better_than_native.csv to get list of complexes to visualize.
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

    import pandas as pd
    
    out_dir = Path(cfg["paths"]["output_folder"])
    base_vis = Path(cfg["paths"]["visuals"])
    rec_dir = Path(cfg["paths"]["receptors_cleaned_folder"])
    crystals_dir = Path(cfg["paths"]["crystals_folder"])
    receptors_dir = Path(cfg["paths"]["receptors_folder"])
    
    # Load better_than_native.csv to get only hits
    better_csv = out_dir / "better_than_native.csv"
    if not better_csv.exists():
        log.warning("better_than_native.csv not found – skipping 2D visualization")
        return
    
    try:
        better_df = pd.read_csv(better_csv)
    except Exception as e:
        log.error("Failed to read better_than_native.csv: %s", e)
        return
    
    if better_df.empty:
        log.info("better_than_native.csv is empty – no complexes to visualize")
        return
    
    log.info("Preparing 2D visualizations for %d complexes from better_than_native.csv…", len(better_df))

    for _, row in better_df.iterrows():
        rec_id = str(row["receptor"])
        lig_id = str(row["ligand"])
        mode = str(row.get("mode", "dock"))
        
        # Map mode to expected tag and result directory
        if mode == "native_redock":
            result_dir = out_dir / "dock_native"
            expected_tag = "native_redock"
        elif mode == "diag":
            result_dir = out_dir / "diagonal"
            expected_tag = "diag"
        else:
            result_dir = out_dir / "matrix"
            expected_tag = "dock"
        
        # Find the PDBQT file
        lig_pdbqt = result_dir / f"{rec_id}__{lig_id}__{expected_tag}.pdbqt"
        if not lig_pdbqt.exists():
            log.warning("PDBQT not found: %s – skipping", lig_pdbqt.name)
            continue

        rec_pdbqt = rec_dir / f"{rec_id}.pdbqt"
        if not rec_pdbqt.exists():
            log.warning("Receptor %s missing – skipping %s", rec_id, lig_pdbqt.name)
            continue

        complex_dir = base_vis / f"{expected_tag}__{rec_id}__{lig_id}"
        complex_dir.mkdir(parents=True, exist_ok=True)

        out_pdb = complex_dir / f"complex__{rec_id}__{lig_id}__{expected_tag}.pdb"
        
        # Try to find original PDB with waters (check both crystals_dir and receptors_dir)
        original_pdb = None
        for candidate_dir in [crystals_dir, receptors_dir]:
            # Try original PDB first
            candidate = candidate_dir / f"{rec_id}.pdb"
            if candidate.exists():
                original_pdb = candidate
                break
            # Try fixed PDB (from PDBFixer)
            candidate = candidate_dir / f"{rec_id}_fixed.pdb"
            if candidate.exists():
                original_pdb = candidate
                break
        
        if not out_pdb.exists():
            try:
                prepare_complex_pdb(rec_pdbqt, lig_pdbqt, complex_dir, original_pdb=original_pdb)
                log.debug("Prepared complex: %s", out_pdb.name)
            except Exception as e:
                log.error("Failed to prepare complex for %s: %s", lig_pdbqt.name, e)
                continue

        plip_dir = complex_dir / "plip"
        plip_xml = plip_dir / "report.xml"
        plip_pdb = complex_dir / f"{rec_id}_{lig_id}_plip.pdb"
        
        # Generate PLIP report if it doesn't exist
        # The improved version automatically includes waters from original PDB
        if not plip_xml.exists():
            # Ensure complex PDB with waters is created
            if not plip_pdb.exists():
                try:
                    prepare_complex_pdb(rec_pdbqt, lig_pdbqt, complex_dir, original_pdb=original_pdb)
                    log.debug("Prepared complex PDB with waters: %s", complex_dir.name)
                except Exception as e:
                    log.warning("Failed to prepare complex PDB for %s: %s", complex_dir.name, e)
                    continue
            
            plip_dir.mkdir(exist_ok=True)
            try:
                # Use PLIP with improved parameters for better interaction detection
                # Note: PLIP automatically detects interactions, but having waters in the structure
                # is crucial for water bridge detection. We've added waters from original PDB above.
                # PLIP will use default thresholds which are well-tuned for most cases.
                plip_cmd = shutil.which("plip") or "plip"
                cmd = [
                    plip_cmd, "-f", str(plip_pdb), "-o", str(plip_dir), "-x"
                ]
                result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                log.debug("Generated PLIP XML for: %s", out_pdb.name)
                if result.stdout:
                    log.debug("PLIP output: %s", result.stdout[:200])
            except subprocess.CalledProcessError as e:
                log.error("PLIP failed for %s: %s", out_pdb.name, e)
                if e.stdout:
                    log.debug("PLIP STDOUT: %s", e.stdout[:500])
                if e.stderr:
                    log.debug("PLIP STDERR: %s", e.stderr[:500])
            except Exception as e:
                log.error("PLIP failed for %s: %s", out_pdb.name, e)


def generate_2d_boards(cfg: dict, log):
    """
    Generate 2D interaction boards for all complexes in better_than_native.csv.

    This function reads better_than_native.csv and generates 2D visualization
    boards for each complex using the visualize_2d module.

    Parameters
    ----------
    cfg : dict
        Pipeline configuration with 'paths'.
    log : Logger
        Logger for progress and errors.

    Returns
    -------
    None
    """
    try:
        from visualize_2d import render_single
    except ImportError:
        log.error("visualize_2d module not found – skipping 2D board generation")
        return

    import pandas as pd
    
    out_dir = Path(cfg["paths"]["output_folder"])
    base_vis = Path(cfg["paths"]["visuals"])
    
    # Load better_than_native.csv
    better_csv = out_dir / "better_than_native.csv"
    if not better_csv.exists():
        log.warning("better_than_native.csv not found – skipping 2D board generation")
        return
    
    try:
        better_df = pd.read_csv(better_csv)
    except Exception as e:
        log.error("Failed to read better_than_native.csv: %s", e)
        return
    
    if better_df.empty:
        log.info("better_than_native.csv is empty – no boards to generate")
        return
    
    log.info("Generating 2D boards for %d complexes from better_than_native.csv…", len(better_df))
    
    boards_dir = base_vis / "2d_boards"
    boards_dir.mkdir(parents=True, exist_ok=True)
    
    success_count = 0
    for idx, row in better_df.iterrows():
        rec_id = str(row["receptor"])
        lig_id = str(row["ligand"])
        mode = str(row.get("mode", "dock"))
        
        # Map mode to expected tag
        if mode == "native_redock":
            expected_tag = "native_redock"
        elif mode == "diag":
            expected_tag = "diag"
        else:
            expected_tag = "dock"
        
        complex_dir = base_vis / f"{expected_tag}__{rec_id}__{lig_id}"
        if not complex_dir.exists():
            log.warning("Complex directory not found: %s – skipping", complex_dir.name)
            continue
        
        # Check if PLIP XML exists
        plip_xml = complex_dir / "plip" / "report.xml"
        if not plip_xml.exists():
            log.warning("PLIP XML not found for %s – skipping", complex_dir.name)
            continue
        
        # Generate 2D board
        out_path = boards_dir / f"{rec_id}__{lig_id}__{expected_tag}.png"
        try:
            render_single(complex_dir, out_path, size=(2200, 1800), fmt="png")
            log.debug("Generated 2D board: %s", out_path.name)
            success_count += 1
        except Exception as e:
            log.error("Failed to generate 2D board for %s: %s", complex_dir.name, e)
            continue
    
    log.info("Generated %d/%d 2D boards successfully", success_count, len(better_df))








