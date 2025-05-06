from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import subprocess, json, time
from utils.validation import files_ok
from tqdm import tqdm
# ── helpers ───────────────────────────────────────────

# ── helpers ───────────────────────────────────────────
def _load_box(receptor_path, cfg):
    box_json = Path(receptor_path).with_suffix('.box.json')
    if box_json.exists():
        with open(box_json) as fh:
            box = json.load(fh)
        return box['center'], box['size']
    # fallback: static box from YAML
    return cfg['docking_params']['center'], cfg['docking_params']['size']


def _single_dock(cfg, rec, lig, is_native):
    """Run one Smina docking. Return (success, log_file, elapsed_s)."""
    outdir = Path(cfg['paths']['output_folder'])
    rec_id, lig_id = Path(rec).stem, Path(lig).stem
    tag = 'native' if is_native else 'dock'
    out_base = outdir / f"{lig_id}__to__{rec_id}__{tag}"

    center, size = _load_box(rec, cfg)
    p = cfg['docking_params']

    cmd = [
        cfg['paths']['smina_path'],
        '--receptor', rec,
        '--ligand',   lig,
        '--out', str(out_base.with_suffix('.pdbqt')),
        '--log', str(out_base.with_suffix('.log')),
        '--exhaustiveness', str(p['exhaustiveness']),
        '--num_modes',      str(p['num_modes']),
        '--cpu',            str(p['cpu']),
        '--center_x', str(center[0]), '--center_y', str(center[1]), '--center_z', str(center[2]),
        '--size_x',   str(size[0]),   '--size_y',   str(size[1]),   '--size_z',   str(size[2]),
    ]
    if is_native:
        cmd += ['--autobox_ligand', lig, '--autobox_add', str(p['autobox_add'])]

    t0 = time.perf_counter()
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elapsed = time.perf_counter() - t0
    ok = files_ok(out_base.parent)    # 3‑file sanity check
    return ok, str(out_base.with_suffix('.log')), elapsed


# ── public API ───────────────────────────────────────
def _build_tasks(recs, ligs, cfg):
    recs = sorted(Path(cfg['paths']['receptors_folder']).glob('*.pdbqt'))
    native_lgs = sorted(Path(cfg['paths']['native_ligands_folder']).glob('*.pdbqt'))
    test_lgs = sorted(Path(cfg['paths']['ligands_folder']).glob('*.pdbqt'))

    # tasks builder
    if cfg['docking_mode'] == 'redock_native':
        tasks = [(cfg, str(r), str(l), True) for r, l in zip(recs, native_lgs)]
    elif cfg['docking_mode'] == 'diagonal':
        tasks = [(cfg, str(r), str(t), False) for r, t in zip(recs, test_lgs)]
    else:  # full_matrix
        tasks = [(cfg, str(r), str(t), False) for r in recs for t in test_lgs]


def run_batch_docking(cfg, log):
    recs = sorted(Path(cfg['paths']['receptors_folder']).glob('*.pdbqt'))
    ligs = sorted(Path(cfg['paths']['ligands_folder']).glob('*.pdbqt'))
    tasks = _build_tasks([str(r) for r in recs], [str(l) for l in ligs], cfg)

    if not tasks:
        log.error("No receptors or ligands found – aborting.")
        return

    log.info(f"{len(tasks)} docking jobs (mode={cfg['docking_mode']})")

    failed, done = [], []
    with ProcessPoolExecutor(max_workers=cfg['runtime']['n_jobs']) as pool:
        for ok, logfile, t in tqdm(pool.map(lambda args: _single_dock(*args), tasks),
                                   total=len(tasks), desc='Docking'):
            log.info(f"{logfile}  –  {'OK' if ok else 'FAIL'}  ({t:.1f}s)")
            (done if ok else failed).append(logfile)

    # ── automatic retry for failed jobs ──
    if failed:
        log.warning(f"Retrying {len(failed)} failed jobs …")
        retry_tasks = [tsk for tsk in tasks
                       if Path(tsk[2]).stem in {Path(f).stem.split('__')[0] for f in failed}]
        with ProcessPoolExecutor(max_workers=cfg['runtime']['n_jobs']) as pool:
            for ok, logfile, t in tqdm(pool.map(lambda args: _single_dock(*args), retry_tasks),
                                       total=len(retry_tasks), desc='Retry'):
                log.info(f"RETRY {logfile}  –  {'OK' if ok else 'FAIL'}  ({t:.1f}s)")