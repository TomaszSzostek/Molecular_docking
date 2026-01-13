
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

from biotite.structure.io import pdb as biotite_pdb
import biotite.structure as biotite_structure

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
        Path to the output processed protein PDB file.

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
    output_filepath = output_dir / f"{input_filepath.stem}.pdb"
    new_input_filepath = output_dir / f"{input_filepath.stem}_old.pdb"
    shutil.move(input_filepath, new_input_filepath)
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

def split_receptor_ligand(
    pdb_path: Path,
    ligand_resname: str,
    receptors_dir: Path,
    ligands_dir: Path,
    overwrite: bool = False
):
    struct = biotite_pdb.get_structure(biotite_pdb.PDBFile.read(pdb_path), include_bonds=True, model=1)

    pdb_id = pdb_path.stem
    rec_path = receptors_dir / f"{pdb_id}.pdb"
    lig_path = ligands_dir / f"{pdb_id}_{ligand_resname}.pdb"

    if not overwrite and rec_path.exists() and lig_path.exists():
        return

    protein_structure = struct[struct.hetero == False]
    ligand_structures = struct[struct.res_name == ligand_resname]

    for ligand_structure in biotite_structure.chain_iter(ligand_structures):
        break

    # Validate ligand = only one residue
    unique_residues = set(ligand_structure.res_id)
    if len(unique_residues) > 1:
        raise RuntimeError(f"❌ Ligand {ligand_resname} in {pdb_path.name} has multiple residues ({len(unique_residues)})")

    ligand_file = biotite_pdb.PDBFile()
    ligand_file.set_structure(ligand_structure)
    ligand_file.write(lig_path)

    protein_file = biotite_pdb.PDBFile()
    protein_file.set_structure(protein_structure)
    protein_file.write(rec_path)

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
