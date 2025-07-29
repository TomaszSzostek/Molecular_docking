"""
rank_vs_native.py

Post‑process docking results to identify test ligands that bind better than
the native co‑crystallized ligand for each receptor, using a user‑defined
affinity margin.
"""

from pathlib import Path
import pandas as pd


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

    # Compute the best (lowest) affinity of native per receptor
    native_scores = (
        native_df
        .groupby("receptor")["min_affinity"]
        .min()
        .rename("native_score")
    )

    # Prepare test ligands by merging with native_score
    hits = (
        df[df["is_native"] == False]
        .merge(native_scores, on="receptor", how="left")
    )

    # Select ligands beating native by at least margin
    hits = hits[hits["min_affinity"] <= hits["native_score"] - margin]

    # Save filtered hits to CSV
    hits_path = out_dir / "better_than_native.csv"
    hits.to_csv(hits_path, index=False)

    log.info(f"Found {len(hits)} ligands better than native.")
    log.info(f"Saved filtered results → {hits_path}")

