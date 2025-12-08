
"""Utility functions for downloading PDB crystal structures and splitting them
into separate receptor/ligand PDB files.

The module is *idempotent*: running it twice will never download or process the
same structure twice. Only receptors/ligands that do **not** already exist on
disk are fetched and generated.

It ensures correct ligand–receptor pairing by extracting the receptor chain
that contains the selected ligand residue (even in symmetric dimers).
"""
from __future__ import annotations
import gzip
import shutil
from pathlib import Path
from typing import Mapping, Sequence

import requests
from Bio.PDB import PDBIO, PDBParser, Select
from prepare_inputs.adt_preparator import ligand_to_pdbqt

try:
    import pdbfixer
    from openmm import app
    PDBFIXER_AVAILABLE = True
except ImportError:
    PDBFIXER_AVAILABLE = False

###############################################################################
# Helpers
###############################################################################

def _guess_biggest_het_ligand(pdb_file: Path) -> str:
    sizes: dict[str, int] = {}
    with open(pdb_file) as fh:
        for line in fh:
            if line.startswith("HETATM"):
                res = line[17:20].strip().upper()
                # Ignore crystallographic waters when guessing ligand
                if res in {"HOH", "WAT"}:
                    continue
                sizes[res] = sizes.get(res, 0) + 1

    if not sizes:
        raise RuntimeError(f"{pdb_file.name} contains no HETATM records – cannot guess ligand.")

    return max(sizes, key=sizes.get)

def preprocess_protein(input_filepath: Path, ph: float = 7.4) -> Path:
    """
    Process a protein file with PDBFixer to prepare it for docking.

    This function uses PDBFixer to add missing residues, atoms, and hydrogens.
    It produces a cleaned and complete protein structure suitable for molecular
    docking simulations.

    Parameters
    ----------
    input_filepath : Path
        Path to the input protein PDB file.
    ph : float
        The pH to use for adding missing hydrogens (default: 7.4).

    Returns
    -------
    Path
        Path to the output processed protein PDB file (suffixed with '_fixed').

    Raises
    ------
    AssertionError
        If input or output file does not exist.
    RuntimeError
        If PDBFixer is not available.
    """
    if not PDBFIXER_AVAILABLE:
        raise RuntimeError("PDBFixer is not installed. Install with: conda install -c conda-forge pdbfixer openmm")
    
    assert input_filepath.exists(), f"PDB file {str(input_filepath)} does not exist!"

    fixer = pdbfixer.PDBFixer(filename=str(input_filepath))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    output_dir = input_filepath.parent
    output_filepath = output_dir / f"{input_filepath.stem}_fixed.pdb"
    with open(output_filepath, 'w') as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    assert output_filepath.exists(), f"Output file {str(output_filepath)} does not exist!"
    return output_filepath


def fetch(pdb_id: str, out_dir: Path, *, overwrite: bool = False) -> Path:
    pdb_id = pdb_id.upper()
    out_dir.mkdir(parents=True, exist_ok=True)

    pdb_path = out_dir / f"{pdb_id}.pdb"
    if pdb_path.exists() and not overwrite:
        return pdb_path

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb.gz"
    gz_path = out_dir / f"{pdb_id}.pdb.gz"

    response = requests.get(url, timeout=30)
    response.raise_for_status()

    gz_path.write_bytes(response.content)
    with gzip.open(gz_path, "rb") as gzf, open(pdb_path, "wb") as out:
        shutil.copyfileobj(gzf, out)
    gz_path.unlink(missing_ok=True)

    return pdb_path

###############################################################################
# Biopython selectors
###############################################################################

class LigandSelect(Select):
    """
    Select only the first full HETATM ligand residue with given resname.
    Also stores which chain the ligand was found in, to extract matching receptor.
    """
    def __init__(self, het_id: str):
        self.het_id = het_id.strip().upper()
        self.target_resid = None
        self.found = 0
        self.selected_chain_id = None

    def accept_residue(self, residue):
        resname = residue.resname.strip().upper()
        hetfield = residue.id[0]
        resid = (residue.get_parent().id, residue.id[1])  # chain ID, residue ID

        if hetfield == " " or resname != self.het_id:
            return False

        if self.target_resid is None:
            self.target_resid = resid
            self.selected_chain_id = resid[0]
            self.found = 1
            return True

        return resid == self.target_resid


class ProteinSelectMatchingChain(Select):
    """Export *only* the chain matching the ligand, excluding water and the ligand."""
    def __init__(self, ligand_resname: str, chain_id: str):
        self.ligand_resname = ligand_resname.strip().upper()
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        hetfield = residue.id[0]
        resname = residue.resname.strip().upper()
        return (
            hetfield == " "
            and resname != self.ligand_resname
            and resname not in {"HOH", "WAT"}
        )

###############################################################################
# Main splitter
###############################################################################

def split_receptor_ligand(
    pdb_path: Path,
    ligand_resname: str,
    receptors_dir: Path,
    ligands_dir: Path,
    overwrite: bool = False
):
    struct = PDBParser(QUIET=True).get_structure('X', str(pdb_path))

    pdb_id = pdb_path.stem
    rec_path = receptors_dir / f"{pdb_id}.pdb"
    lig_path = ligands_dir / f"{pdb_id}_{ligand_resname}.pdb"

    if not overwrite and rec_path.exists() and lig_path.exists():
        return

    io = PDBIO()

    # Save ligand (first matching resname only)
    io.set_structure(struct)
    lig_sel = LigandSelect(ligand_resname)
    io.save(str(lig_path), lig_sel)

    if lig_sel.found == 0:
        rec_path.unlink(missing_ok=True)
        lig_path.unlink(missing_ok=True)
        raise RuntimeError(f"❌ Ligand '{ligand_resname}' not found in {pdb_path.name}")

    # Save receptor (matching ligand's chain)
    io.set_structure(struct)
    io.save(str(rec_path), ProteinSelectMatchingChain(ligand_resname, lig_sel.selected_chain_id))

    # Validate ligand = only one residue
    residues = set()
    with open(lig_path) as f:
        for line in f:
            if line.startswith("HETATM"):
                resid = (line[21], line[22:26])
                residues.add(resid)
    if len(residues) > 1:
        rec_path.unlink(missing_ok=True)
        lig_path.unlink(missing_ok=True)
        raise RuntimeError(f"❌ Ligand {ligand_resname} in {pdb_path.name} has multiple residues ({len(residues)})")

    for extra in receptors_dir.glob(f"{pdb_id}_*.pdb"):
        extra.unlink(missing_ok=True)

###############################################################################
# Batch API used by the pipeline
###############################################################################

def fetch_and_split_batch(cfg: Mapping, log, *, overwrite: bool = False) -> None:
    pdb_ids: Sequence[str] = cfg["fetch_crystals"]["pdb_ids"]
    res_map: Mapping[str, str] = cfg["fetch_crystals"].get("ligand_resnames", {})

    raw_dir = Path(cfg["paths"]["crystals_folder"]).expanduser()
    receptors_dir = Path(cfg["paths"]["receptors_folder"]).expanduser()
    ligands_dir = Path(cfg["paths"]["native_ligands_folder"]).expanduser()

    done_recs = {p.stem.upper() for p in receptors_dir.glob("*.pdb")}
    done_ligs = {p.stem.split("_")[0].upper() for p in ligands_dir.glob("*.pdb")}

    for pdb_id in pdb_ids:
        pdb_id_upper = pdb_id.upper()

        if (
            pdb_id_upper in done_recs
            and pdb_id_upper in done_ligs
            and not overwrite
        ):
            log.info("✓ %s already processed – skipping", pdb_id)
            continue

        try:
            log.info("Downloading %s…", pdb_id)
            pdb_path = fetch(pdb_id_upper, raw_dir, overwrite=overwrite)

            # Preprocess with PDBFixer if available
            if PDBFIXER_AVAILABLE:
                log.info("Preprocessing %s with PDBFixer…", pdb_id)
                try:
                    pdb_path = preprocess_protein(pdb_path)
                    log.info("✓ PDBFixer complete → %s", pdb_path.name)
                except Exception as exc:
                    log.warning("PDBFixer failed for %s: %s, using original", pdb_id, exc)

            ligand_resname = res_map.get(pdb_id)
            if not ligand_resname:
                ligand_resname = _guess_biggest_het_ligand(pdb_path)
                log.warning("No ligand_resname for %s, guessed '%s'", pdb_id, ligand_resname)

            log.info("Splitting %s → receptor + %s", pdb_id, ligand_resname)
            split_receptor_ligand(
                pdb_path,
                ligand_resname,
                receptors_dir,
                ligands_dir,
                overwrite=overwrite,
            )
        except Exception as exc:
            log.error("Failed processing %s: %s", pdb_id, exc)
            continue
