import logging, os, sys
from pathlib import Path


def ensure_dirs(cfg):
    for k in ['receptors_folder','ligands_folder','output_folder','temp_dir','crystals_folder']:
        Path(cfg['paths'][k]).mkdir(parents=True, exist_ok=True)


def setup_logging(cfg):
    log_path = Path(cfg['paths']['output_folder']) / 'run.log'
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s | %(levelname)s | %(message)s',
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger('dockpipe')