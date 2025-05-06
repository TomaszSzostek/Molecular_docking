#!/usr/bin/env python
import argparse, yaml, sys
from utils.file_utils import setup_logging, ensure_dirs
from dock.docking import run_batch_docking
from analyze.results_extractor import consolidate_logs
from analyze.postprocess import rank_vs_native


def load_config(cfg_path):
    with open(cfg_path) as fh:
        return yaml.safe_load(fh)


def cli():
    p = argparse.ArgumentParser()
    p.add_argument('--config', default='config.yaml', help='YAML config file')
    return p.parse_args()


def main():
    args = cli()
    cfg = load_config(args.config)
    ensure_dirs(cfg)  # <-- creates data/* paths if absent
    log = setup_logging(cfg)

    log.info('>>> START DOCKING PIPELINE')
    run_batch_docking(cfg, log)

    log.info('>>> EXTRACT & MERGE LOGS')
    consolidate_logs(cfg, log)

    log.info('>>> POST‑PROCESS & RANK')
    rank_vs_native(cfg, log)

    log.info('Pipeline finished ✓')


if __name__ == '__main__':
    sys.exit(main())