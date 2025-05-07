import argparse, yaml, sys
from pathlib import Path
from utils.file_utils import setup_logging, ensure_dirs
from utils.validation import files_ok
from prepare.adt_preparator import receptor_to_pdbqt, ligand_to_pdbqt
from prepare.fetch_crystals import fetch_and_split_batch
from prepare.smiles_loader import csv_to_pdbqt
from dock.docking import run_batch_docking
from analyze.results_extractor import consolidate_logs
from analyze.postprocess import rank_vs_native
from analyze.RMSD import run_rmsd_and_plot


def load_config(p):
    return yaml.safe_load(open(p))

def cli():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default="config.yaml")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()

def prepare_inputs(cfg, log):
    rec_dir = Path(cfg["paths"]["receptors_folder"])
    rec_clean = Path(cfg["paths"]["receptors_cleaned_folder"])
    nat_dir = Path(cfg["paths"]["native_ligands_folder"])
    test_dir = Path(cfg["paths"]["ligands_folder"])

    # Skip if cleaned PDBQT receptors already exist
    if not list(rec_clean.glob("*.pdbqt")):
        for pdb in rec_dir.glob("*.pdb"):
            receptor_to_pdbqt(pdb, rec_clean)

    # Native ligands (sanitization OFF)
    if not list(nat_dir.glob("*.pdbqt")):
        for pdb in nat_dir.glob("*.pdb"):
            ligand_to_pdbqt(pdb, nat_dir, sanitize=False)

    # Test ligands (sanitization ON)
    if not list(test_dir.glob("*.pdbqt")):
        for pdb in test_dir.glob("*.pdb"):
            ligand_to_pdbqt(pdb, test_dir, sanitize=True)

    log.info("PDB → PDBQT conversion finished")


def main():
    args = cli()
    cfg = load_config(args.config)
    ensure_dirs(cfg)
    log = setup_logging(cfg)

    if cfg.get("fetch_crystals", {}).get("enabled"):
        fetch_and_split_batch(cfg, log)

    if cfg.get("ligand_csv", {}).get("enabled"):
        smi_cfg = cfg["ligand_csv"]
        csv_to_pdbqt(Path(smi_cfg["file"]),
                     smi_cfg["smiles_column"],
                     smi_cfg["id_column"],
                     Path(cfg["paths"]["ligands_folder"]))

    if args.dry_run:
        from dock.docking import _build_tasks
        rec = list(Path(cfg["paths"]["receptors_cleaned_folder"]).glob("*.pdbqt"))
        lig = list(Path(cfg["paths"]["ligands_folder"]).glob("*.pdbqt"))
        print(f"[DRY‑RUN] would dock {len(rec)} receptors × {len(lig)} ligands")
        sys.exit(0)

    log.info(">>> PREPARE INPUTS")
    prepare_inputs(cfg, log)

    log.info(">>> DOCK")
    out_dir = Path(cfg["paths"]["output_folder"])
    existing_docked = list(out_dir.glob("*redock.pdbqt"))
    if cfg["docking_mode"] == "redock_native" and existing_docked:
        log.info("Docking already performed. Skipping.")
    else:
        run_batch_docking(cfg, log)

    log.info(">>> MERGE LOGS")
    consolidate_logs(cfg, log)

    log.info(">>> POST‑PROCESS")
    rank_vs_native(cfg, log)

    if cfg["docking_mode"] == "redock_native":
        log.info(">>> RMSD & VISUALIZATION")
        run_rmsd_and_plot(cfg, log)

    out = Path(cfg["paths"]["output_folder"])
    ok = any("native_redock" in p.name for p in out.glob("*.pdbqt"))
    log.info("Pipeline finished – validation %s", "OK" if ok else "FAIL")


if __name__ == "__main__":
    sys.exit(main())
