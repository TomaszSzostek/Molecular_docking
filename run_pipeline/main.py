"""
Main pipeline script for molecular docking using Smina.

This script performs the following steps:
1. PREPARE INPUTS – Fetch PDBs and convert receptors/ligands to PDBQT.
2. DOCK – Run docking simulations using Smina.
3. MERGE LOGS – Consolidate docking output logs.
4. POST‑PROCESSING – Rank results and compare to native ligands.
5. RMSD & OVERLAYS – Calculate RMSD and generate 2D overlays (only in redock mode).
6. FINAL CHECK – Validate success of docking pipeline.
7. VISUALIZATION – Preparing files for visualization in jupyter notebooks.
8. (Optional) Convert SMILES from CSV to ligands.
9. DONE – Complete.

This script reads from a YAML configuration file and supports dry-run mode.
"""

import argparse
import yaml
import sys
from pathlib import Path

from utils.file_utils import setup_logging, ensure_dirs
from utils.validation import files_ok
from prepare_inputs.adt_preparator import receptor_to_pdbqt, ligand_to_pdbqt
from prepare_inputs.fetch_crystals import fetch_and_split_batch
from prepare_inputs.smiles_loader import csv_to_pdbqt
from dock_smina.docking import run_batch_docking
from analyze.results_extractor import consolidate_logs
from analyze.postprocess import rank_vs_native
from analyze.RMSD import run_rmsd_and_plot
from analyze.files_for_visualization import generate_files


def load_config(path):
    """Load YAML configuration file."""
    return yaml.safe_load(open(path))


def cli():
    """Command-line interface."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def prepare_inputs(cfg, log):
    """
    Step 1/9: Prepare receptor and ligand inputs.
    Converts .pdb to .pdbqt and fetches PDB structures if needed.
    """
    log.info("[1/9] PREPARE INPUTS")

    rec_dir = Path(cfg["paths"]["receptors_folder"])
    rec_clean = Path(cfg["paths"]["receptors_cleaned_folder"])
    nat_dir = Path(cfg["paths"]["native_ligands_folder"])
    test_dir = Path(cfg["paths"]["ligands_folder"])

    if not list(rec_dir.glob("*.pdb")):
        log.info("→ Fetching crystal structures from PDB...")
        fetch_and_split_batch(cfg, log)
    else:
        log.info("→ Crystal structures already exist – skipping fetch.")

    if not list(rec_clean.glob("*.pdbqt")):
        log.info("→ Converting receptors to PDBQT...")
        for pdb in rec_dir.glob("*.pdb"):
            receptor_to_pdbqt(pdb, rec_clean)
    else:
        log.info("→ Receptor PDBQT files already exist – skipping conversion.")

    if not list(nat_dir.glob("*.pdbqt")):
        log.info("→ Converting native ligands to PDBQT...")
        for pdb in nat_dir.glob("*.pdb"):
            ligand_to_pdbqt(pdb, nat_dir, sanitize=False)
    else:
        log.info("→ Native ligand PDBQT files already exist – skipping conversion.")

    if not list(test_dir.glob("*.pdbqt")):
        log.info("→ Converting test ligands to PDBQT...")
        for pdb in test_dir.glob("*.pdb"):
            ligand_to_pdbqt(pdb, test_dir, sanitize=True)
    else:
        log.info("→ Test ligand PDBQT files already exist – skipping conversion.")

    log.info("→ PDB → PDBQT conversion finished")


def main():
    args = cli()
    cfg = load_config(args.config)

    for folder in cfg['paths'].values():
        if isinstance(folder, str) and not folder.startswith('/'):
            Path(folder).mkdir(parents=True, exist_ok=True)

    # Determine docking mode
    mode = cfg.get("docking_mode", "matrix")

    if mode == "redock_native":
        cfg["paths"]["active_ligands_folder"] = cfg["paths"]["native_ligands_folder"]
        cfg["docking_params"]["num_modes"] = 1
    elif mode in {"diagonal", "matrix", "full_matrix"}:
        cfg["paths"]["active_ligands_folder"] = cfg["paths"]["ligands_folder"]
    else:
        raise ValueError(f"Unsupported docking_mode: {mode}")

    ensure_dirs(cfg)
    log = setup_logging(cfg)

    # Optional SMILES to PDBQT conversion
    if cfg.get("ligand_csv", {}).get("enabled"):
        smi_cfg = cfg["ligand_csv"]
        csv_to_pdbqt(Path(smi_cfg["file"]),
                     smi_cfg["smiles_column"],
                     smi_cfg["id_column"],
                     Path(cfg["paths"]["ligands_folder"]))

    if args.dry_run:
        rec = list(Path(cfg["paths"]["receptors_cleaned_folder"]).glob("*.pdbqt"))
        lig = list(Path(cfg["paths"]["active_ligands_folder"]).glob("*.pdbqt"))
        print(f"[DRY‑RUN] would dock_smina {len(rec)} receptors × {len(lig)} ligands")
        sys.exit(0)

    prepare_inputs(cfg, log)

    # Step 2/9: Docking
    log.info("[2/9] DOCK")
    out_root = Path(cfg["paths"]["output_folder"]).resolve()
    mode = cfg.get("docking_mode", "matrix")
    if mode == "redock_native":
        out_dir = out_root / "dock_native"
    elif mode == "diagonal":
        out_dir = out_root / "diagonal"
    elif mode in {"matrix", "full_matrix"}:
        out_dir = out_root / "matrix"
    else:
        raise ValueError(f"Unsupported docking_mode: {mode}")
    already_docked = any(out_dir.glob("*__native_redock__*.pdbqt"))
    if cfg["docking_mode"] == "redock_native" and already_docked:
        log.info("→ Docking already performed – skipping.")
    else:
        run_batch_docking(cfg, log)

    # Step 3/9: Merge docking logs
    log.info("[3/9] MERGE LOGS")
    consolidate_logs(cfg, log)

    # Step 4/9: Post-processing
    log.info("[4/9] POST‑PROCESSING")
    rank_vs_native(cfg, log)

    # Step 5/9: RMSD & overlays (redocking only)
    if cfg["docking_mode"] == "redock_native":
        log.info("[5/9] RMSD & OVERLAYS")
        run_rmsd_and_plot(cfg, log)

    # Step 6/9: Final validation
    log.info("[6/9] FINAL CHECK")
    result_files = list(out_dir.glob("*.pdbqt"))
    successes = sum(1 for f in result_files if f.stat().st_size > 0)
    total = len(result_files)
    ok = successes > 0
    log.info(f"→ Found {successes}/{total} valid docking results.")
    log.info("→ Pipeline finished – validation %s", "OK" if ok else "FAIL")

    # Step 7/9: Visual summary

    log.info("[7/9] Preparing files for visualization")
    generate_files(cfg, log)

    log.info("[9/9] DONE. Full docking workflow completed.")



if __name__ == "__main__":
    main()
