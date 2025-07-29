"""
prepare_inputs/smiles_loader.py

Convert a CSV file of SMILES strings into 3D PDBQT ligand files.

This module provides:
  - Robust CSV reading with BOM handling and auto delimiter detection.
  - Unique ID assignment to avoid filename collisions.
  - Main function to generate .pdbqt files using RDKit and AutoDockTools.

Public API:
    csv_to_pdbqt(csv_path: Path, smiles_col: str, id_col: str, out_dir: Path) -> None
"""

from __future__ import annotations

from pathlib import Path
import logging
import pandas as pd
from rdkit import Chem
from typing import Iterable

from prepare_inputs.adt_preparator import ligand_to_pdbqt

log = logging.getLogger(__name__)

# ───────────────────────── helpers ─────────────────────────

def _read_csv_safely(csv_path: Path, id_col: str, smiles_col: str) -> pd.DataFrame:
    """
    Safely read a CSV file and return a DataFrame containing only valid ID/SMILES rows.

    The function:
      1. Reads with 'utf-8-sig' to consume any BOM.
      2. Attempts auto delimiter detection, then falls back to ';' and tab.
      3. Strips whitespace and BOM characters from column headers.
      4. Validates presence of id_col and smiles_col (case-insensitive).
      5. Drops rows with empty or missing ID or SMILES values.

    Parameters
    ----------
    csv_path : Path
        Path to the input CSV file.
    id_col : str
        Name of the column for ligand identifiers.
    smiles_col : str
        Name of the column for SMILES strings.

    Returns
    -------
    pd.DataFrame
        DataFrame filtered to valid (ID, SMILES) pairs.

    Raises
    ------
    ValueError
        If the file cannot be read, required columns are missing, or no valid rows remain.
    """
    enc = "utf-8-sig"  # use encoding that strips BOM if present
    tried: list[tuple[str | None, str]] = [
        (None, "auto"),  # let pandas guess
        (";", ";"),      # semicolon-separated
        ("\t", "\\t"),   # tab-separated
    ]
    last_exc: Exception | None = None
    df = None

    # Attempt to read with different separators
    for sep, label in tried:
        try:
            df = pd.read_csv(csv_path, sep=sep, engine="python", encoding=enc)
            log.debug("Read CSV (%s sep). Columns: %s", label, df.columns.tolist())
            break
        except Exception as e:
            last_exc = e

    if df is None:
        raise ValueError(f"Could not read CSV {csv_path}: {last_exc}")

    # Normalize header names: strip whitespace and BOM
    df.rename(columns=lambda c: str(c).strip().lstrip("\ufeff"), inplace=True)

    # If required columns missing, try case-insensitive match
    cols_lower = {c.lower(): c for c in df.columns}
    if id_col not in df.columns or smiles_col not in df.columns:
        want = {id_col.lower(): id_col, smiles_col.lower(): smiles_col}
        for low, orig in want.items():
            if low in cols_lower:
                df.rename(columns={cols_lower[low]: orig}, inplace=True)

    # Final check for required columns
    if id_col not in df.columns or smiles_col not in df.columns:
        raise ValueError(
            f"CSV must contain columns '{id_col}' and '{smiles_col}', "
            f"but found: {df.columns.tolist()}"
        )

    # Clean and drop rows with empty or missing values
    df[id_col] = df[id_col].astype(str).str.strip()
    df[smiles_col] = df[smiles_col].astype(str).str.strip()
    df = df[(df[id_col] != "") & (df[smiles_col] != "")]
    df.dropna(subset=[id_col, smiles_col], inplace=True)

    if df.empty:
        raise ValueError("CSV has no valid (ID, SMILES) rows after cleaning.")

    return df


def _uniqueify_ids(ids: Iterable[str], existing: set[str]) -> list[str]:
    """
    Ensure each ID is unique by appending numeric suffixes to duplicates.

    Parameters
    ----------
    ids : Iterable[str]
        Original sequence of IDs from the CSV.
    existing : set[str]
        IDs (filenames) that already exist in the output directory.

    Returns
    -------
    list[str]
        New list of IDs guaranteed not to collide, e.g. ['lig', 'lig_1', ...].
    """
    out = []
    seen = set(existing)
    for raw in ids:
        base = raw
        k = 0
        new_id = base
        # Append suffix until a unique ID is found
        while new_id in seen:
            k += 1
            new_id = f"{base}_{k}"
        seen.add(new_id)
        out.append(new_id)
    return out

# ───────────────────────── main ─────────────────────────

def csv_to_pdbqt(csv_path: Path,
                 smiles_col: str,
                 id_col: str,
                 out_dir: Path) -> None:
    """
    Read a CSV of SMILES strings and generate PDBQT ligand files.

    Workflow:
      1. Read and clean CSV via _read_csv_safely.
      2. Create output directory if it does not exist.
      3. Assign unique ligand IDs to avoid overwriting existing files.
      4. For each row:
         a. Parse SMILES into RDKit Mol; skip invalid SMILES.
         b. Add hydrogens and write temporary PDB.
         c. Convert temporary PDB to PDBQT via ligand_to_pdbqt.
         d. Remove temporary PDB file.
      5. Log counts of successful and skipped conversions.

    Parameters
    ----------
    csv_path : Path
        Path to the input CSV file.
    smiles_col : str
        Name of the SMILES column.
    id_col : str
        Name of the ligand ID column.
    out_dir : Path
        Directory in which to write <id>.pdbqt files.

    Returns
    -------
    None
    """
    # Load and clean the input CSV
    df = _read_csv_safely(csv_path, id_col=id_col, smiles_col=smiles_col)

    # Ensure the output directory exists
    out_dir.mkdir(parents=True, exist_ok=True)

    # Determine existing ligand filenames to avoid collisions
    existing = {p.stem for p in out_dir.glob("*.pdbqt")}
    df[id_col] = _uniqueify_ids(df[id_col].tolist(), existing)

    ok, bad = 0, 0
    for _, row in df.iterrows():
        ligand_id = row[id_col]
        smiles    = row[smiles_col]

        # Parse SMILES string into molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            log.warning("Invalid SMILES for id=%s: %s – skipping", ligand_id, smiles)
            bad += 1
            continue

        # Add explicit hydrogens for 3D conversion
        mol = Chem.AddHs(mol)

        tmp_pdb = out_dir / f"{ligand_id}.pdb"
        try:
            # Write temporary PDB file
            Chem.MolToPDBFile(mol, str(tmp_pdb))
        except Exception as e:
            log.error("Failed to write PDB for %s: %s", ligand_id, e)
            bad += 1
            if tmp_pdb.exists():
                tmp_pdb.unlink()
            continue

        # Convert the PDB to PDBQT using AutoDockTools
        try:
            ligand_to_pdbqt(tmp_pdb, out_dir, sanitize=True)
            ok += 1
        except Exception as e:
            log.error("Failed to convert %s to PDBQT: %s", tmp_pdb.name, e)
            bad += 1
        finally:
            # Clean up the temporary PDB file
            if tmp_pdb.exists():
                tmp_pdb.unlink()

    log.info("SMILES→PDBQT finished: %d ok, %d skipped.", ok, bad)

