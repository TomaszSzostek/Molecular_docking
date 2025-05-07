import re
import pandas as pd
from pathlib import Path

aff_pattern = re.compile(r"^\s*\d+\s+(-?\d+\.\d+)\s+([\d\.]+)\s+([\d\.]+)")

def parse_smina_log(log_path: Path) -> dict:
    lines = log_path.read_text().splitlines()
    affinities = []

    for line in lines:
        m = aff_pattern.match(line)
        if m:
            affinity = float(m.group(1))
            affinities.append(affinity)

    parts = log_path.stem.split("__")
    receptor = parts[0]
    ligand = parts[1] if len(parts) > 1 else "UNKNOWN"
    mode = parts[2] if len(parts) > 2 else "unknown"

    return {
        "receptor": receptor,
        "ligand": ligand,
        "mode": mode,
        "min_affinity": min(affinities) if affinities else None,
        "max_affinity": max(affinities) if affinities else None,
        "log_path": str(log_path),
        "pose_path": str(log_path.with_suffix(".pdbqt")),
        "is_native": mode == "native_redock"
    }


def consolidate_logs(cfg, log):
    out = Path(cfg['paths']['output_folder'])
    logs = [f for f in out.glob("*.log") if "__" in f.name and f.name.endswith(".log")]
    log.info(f'Parsing {len(logs)} logs …')
    recs = [parse_smina_log(f) for f in logs]
    df = pd.DataFrame(recs)
    df.dropna(subset=["min_affinity"], inplace=True)
    df.to_csv(out / "results.csv", index=False)


