"""
rank_vs_native.py

Post‑process docking results to identify test ligands that bind better than
the native co‑crystallized ligand for each receptor, using a user‑defined
affinity margin.
"""

from pathlib import Path
from typing import Optional, Set

import numpy as np
import pandas as pd

from analyze import RMSD

RMSD_THRESHOLD = 2.0


def rank_vs_native(cfg: dict, log) -> None:
    """
    Filter and save ligands whose best docking score beats the native ligand.

    This function reads the consolidated `results.csv`, verifies that it
    contains native redocking data, computes the minimum affinity of the
    native ligand per receptor, and selects test ligands whose affinity is
    at least `margin` kcal/mol better. The hits are written to
    `better_than_native.csv`.

    Parameters
    ----------
    cfg : dict
        Pipeline configuration dictionary. Expects:
          - paths.output_folder: base output directory
          - postprocess.native_margin: affinity margin threshold (float)
    log : logging.Logger
        Logger for info/warning/error messages.

    Returns
    -------
    None
    """
    out_dir = Path(cfg["paths"]["output_folder"])
    results_file = out_dir / "results.csv"

    # Ensure docking results exist
    if not results_file.exists():
        log.error("results.csv not found. Run docking first.")
        return

    # Load full results table
    df = pd.read_csv(results_file)

    # Affinity margin (kcal/mol) by which test ligands must improve
    margin = cfg["postprocess"].get("native_margin", 0.0)

    # Skip if not in native redock mode
    first_mode = df["mode"].iloc[0]
    if not (first_mode.startswith("native_redock") or "redock" in first_mode):
        log.info("Skipping postprocess – not native redocking mode.")
        return

    # Require presence of 'is_native' flag to distinguish ligands
    if "is_native" not in df.columns:
        log.error("Missing 'is_native' column in results.csv. Cannot rank hits.")
        return

    # Extract native ligand entries
    native_df = df[df["is_native"] == True]
    if native_df.empty:
        log.warning("No native ligands found in results.csv.")
        return

    # Compute the best (lowest) and worst (highest) affinity of native per receptor
    native_scores = (
        native_df
        .groupby("receptor")
        .agg({
            "min_affinity": "min",  # best pose
            "max_affinity": "max"    # worst pose
        })
        .rename(columns={
            "min_affinity": "native_score",
            "max_affinity": "native_max_affinity"
        })
    )

    # Prepare test ligands by merging with native_score
    hits = (
        df[df["is_native"] == False]
        .merge(native_scores, on="receptor", how="left")
    )
    rmsd_map = _native_rmsd_by_receptor(cfg, log)
    if rmsd_map is not None:
        hits = hits.merge(rmsd_map, on="receptor", how="left")
        hits = hits[hits["native_rmsd"].notna()]
        hits = hits[hits["native_rmsd"] <= RMSD_THRESHOLD]
    else:
        log.warning("RMSD summary missing – cannot filter receptors by <=2Å.")

    # Select ligands beating native by at least margin
    hits = hits[hits["min_affinity"] <= hits["native_score"] - margin]

    # Compute pose stability RMSD for each hit
    hits["pose_stability_rmsd"] = hits.apply(
        lambda row: _pose_stability_rmsd(row, out_dir, log), axis=1
    )

    # Drop internal-only columns
    hits = hits.drop(columns=["is_native"], errors="ignore")

    # Save filtered hits to CSV
    hits_path = out_dir / "better_than_native.csv"
    hits.to_csv(hits_path, index=False)

    log.info(f"Found {len(hits)} ligands better than native with stable natives.")
    log.info(f"Saved filtered results → {hits_path}")


def _native_rmsd_by_receptor(cfg: dict, log) -> Optional[pd.DataFrame]:
    """
    Load native redock RMSD summary and map receptor -> RMSD value.
    """
    # RMSD summary is stored in the main results folder (sibling of visuals)
    summary_fp = Path(cfg["paths"]["output_folder"]) / "rmsd_summary.csv"
    if not summary_fp.exists():
        log.warning("rmsd_summary.csv not found at %s", summary_fp)
        return None
    df = pd.read_csv(summary_fp)
    
    rmsd_col = "RMSD_best" if "RMSD_best" in df.columns else "RMSD"
    if "Complex" not in df.columns or rmsd_col not in df.columns:
        log.warning("rmsd_summary.csv missing required columns.")
        return None
    def _extract_receptor(label: str) -> Optional[str]:
        """
        Extract receptor ID from Complex label '<rec>__<lig>__native_redock'.
        """
        parts = str(label).split("__")
        return parts[0] if len(parts) >= 3 else None
    df["receptor"] = df["Complex"].apply(_extract_receptor)
    df = df.dropna(subset=["receptor"])
    if df.empty:
        return None
    mapped = (
        df.groupby("receptor")[rmsd_col]
        .min()
        .reset_index()
        .rename(columns={rmsd_col: "native_rmsd"})
    )
    return mapped


def _pose_stability_rmsd(row, out_dir: Path, log) -> Optional[float]:
    """
    Calculate RMSD spread of docking poses for a ligand (lower = more stable).
    """
    pdbqt_fp = _find_pose_file(out_dir, row["receptor"], row["ligand"], row["mode"])
    if not pdbqt_fp:
        log.debug("Pose file not found for %s__%s__%s", row["receptor"], row["ligand"], row["mode"])
        return None
    try:
        models = RMSD.read_pdbqt_models(pdbqt_fp)
    except Exception as exc:
        log.debug("Failed to read models from %s: %s", pdbqt_fp, exc)
        return None
    if len(models) <= 1:
        return 0.0
    ref_names, _, ref_xyz = models[0]
    rms_values = []
    for names, _, xyz in models[1:]:
        amap = RMSD._map_by_name(ref_names, names)
        if len(amap) < RMSD.MIN_COMMON_ATOMS:
            continue
        p_idx, r_idx = zip(*amap)
        rmsd, _ = RMSD._kabsch(xyz[list(p_idx)], ref_xyz[list(r_idx)])
        rms_values.append(rmsd)
    if not rms_values:
        return 0.0
    return float(np.sqrt(np.mean(np.square(rms_values))))


def _find_pose_file(out_dir: Path, receptor: str, ligand: str, mode: str) -> Optional[Path]:
    """
    Locate the PDBQT file storing all poses for (receptor, ligand, mode).
    """
    pattern = f"{receptor}__{ligand}__{mode}.pdbqt"
    matches = list(out_dir.rglob(pattern))
    return matches[0] if matches else None