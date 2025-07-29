"""
results_extractor.py

Parse Smina log files to extract docking affinities and consolidate them
into a single CSV results table.
"""

import re
from pathlib import Path
import pandas as pd

# Regular expression to match affinity lines in Smina log:
#   e.g. "    1   -7.23   0.000   0.000"
AFF_PATTERN = re.compile(r"^\s*\d+\s+(-?\d+\.\d+)\s+[\d\.]+\s+[\d\.]+")


def parse_smina_log(log_path: Path) -> dict:
    """
    Extract docking affinity statistics from a single Smina log file.

    This reads all lines matching AFF_PATTERN, converts them to floats,
    and determines the minimum and maximum affinities. It also infers
    receptor, ligand, mode, and whether this is a native redock run
    from the filename.

    Parameters
    ----------
    log_path : pathlib.Path
        Path to a .log file named '<receptor>__<ligand>__<mode>.log'.

    Returns
    -------
    dict
        A record containing:
          - receptor (str): receptor identifier
          - ligand (str): ligand identifier (or "UNKNOWN")
          - mode (str): docking mode tag
          - is_native (bool): True if mode contains "native"
          - min_affinity (float or None): best (lowest) affinity, or None
          - max_affinity (float or None): worst (highest) affinity, or None
    """
    # Read all lines and find those matching the affinity pattern
    lines = log_path.read_text().splitlines()
    affinities = [
        float(m.group(1))
        for line in lines
        if (m := AFF_PATTERN.match(line))
    ]

    # Parse filename parts: receptor, ligand, mode
    parts = log_path.stem.split("__")
    receptor = parts[0]
    ligand = parts[1] if len(parts) > 1 else "UNKNOWN"
    mode = parts[2] if len(parts) > 2 else "unknown"

    # Determine if this is a native redock run
    is_native = "native" in mode.lower()

    # Assign min/max affinities
    if affinities:
        if is_native:
            # In native redock, only the first affinity is relevant
            min_aff = max_aff = affinities[0]
        else:
            min_aff, max_aff = min(affinities), max(affinities)
    else:
        # No affinities parsed
        min_aff = max_aff = None

    return {
        "receptor": receptor,
        "ligand": ligand,
        "mode": mode,
        "is_native": is_native,
        "min_affinity": min_aff,
        "max_affinity": max_aff,
    }


def consolidate_logs(cfg: dict, log) -> None:
    """
    Locate all Smina .log files under the results folder, parse them,
    and write a consolidated `results.csv`.

    This function:
      1. Recursively finds '*.log' files containing '__' in their names.
      2. Parses each file to extract affinity stats.
      3. Builds a DataFrame, computes an affinity range string.
      4. Drops entries with no affinity data.
      5. Saves to 'results.csv' in the output folder.

    Parameters
    ----------
    cfg : dict
        Pipeline configuration dictionary. Expects:
          - paths.output_folder: base output directory (str)
    log : logging.Logger
        Logger for informational and error messages.

    Returns
    -------
    None
    """
    out_dir = Path(cfg["paths"]["output_folder"])

    # Gather all .log files produced by Smina
    logs = list(out_dir.rglob("*.log"))
    logs = [f for f in logs if "__" in f.name]
    log.info("Parsing %d log files …", len(logs))

    if not logs:
        log.warning("No .log files found – skipping consolidation.")
        return

    # Parse each log into a dict record
    records = [parse_smina_log(f) for f in logs]
    df = pd.DataFrame(records)

    # Create a human-readable affinity range column
    df["affinity_range"] = df.apply(
        lambda row: f"{row['min_affinity']} | {row['max_affinity']}"
        if pd.notnull(row["min_affinity"]) and pd.notnull(row["max_affinity"])
        else "NA",
        axis=1
    )

    # Remove any rows lacking a min_affinity value
    if "min_affinity" in df.columns:
        df = df.dropna(subset=["min_affinity"], how="all")

    # Write the consolidated results to CSV
    out_file = out_dir / "results.csv"
    df.to_csv(out_file, index=False)
    log.info("Saved consolidated results → %s", out_file)


