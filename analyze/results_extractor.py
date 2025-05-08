import re
from pathlib import Path
import pandas as pd

# 1) szukamy linii w formacie „  1   -7.23   0.000   0.000”
aff_pattern = re.compile(r"^\s*\d+\s+(-?\d+\.\d+)\s+[\d\.]+\s+[\d\.]+")

def parse_smina_log(log_path: Path) -> dict:
    lines = log_path.read_text().splitlines()
    affinities = [
        float(m.group(1))
        for line in lines
        if (m := aff_pattern.match(line))
    ]

    parts = log_path.stem.split("__")
    receptor = parts[0]
    ligand   = parts[1] if len(parts) > 1 else "UNKNOWN"
    mode     = parts[2] if len(parts) > 2 else "unknown"

    is_native = mode in {"native_redock", "dock_native"}

    if affinities:
        if is_native:
            single = affinities[0]
            min_aff = max_aff = single
        else:
            min_aff, max_aff = min(affinities), max(affinities)
    else:
        # brak dopasowań w logu
        min_aff = max_aff = None

    return {
        "receptor":      receptor,
        "ligand":        ligand,
        "mode":          mode,
        "is_native":     is_native,
        "affinity":      min_aff,   # alias – dla wygody
        "min_affinity":  min_aff,
        "max_affinity":  max_aff,
        "log_path":      str(log_path),
        "pose_path":     str(log_path.with_suffix(".pdbqt")),
    }


def consolidate_logs(cfg, log):

    out_dir = Path(cfg["paths"]["output_folder"])
    logs = [f for f in out_dir.glob("*.log") if "__" in f.name]
    log.info("Parsing %d logs …", len(logs))

    if not logs:
        log.warning("No *.log files found – skipping consolidation")
        return

    # Konwersja logów → słowniki → DataFrame
    records = [parse_smina_log(f) for f in logs]
    df = pd.DataFrame(records)

    # Zostawiamy tylko te wiersze, w których jest sensowna wartość
    # (parser zawsze tworzy min_affinity; w trybie native = affinity = min = max)
    if "min_affinity" in df.columns:
        df = df.dropna(subset=["min_affinity"], how="all")

    out_file = out_dir / "results.csv"
    df.to_csv(out_file, index=False)
    log.info("Saved consolidated results → %s", out_file)


