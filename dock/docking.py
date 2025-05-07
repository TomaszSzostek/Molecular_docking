import subprocess
from multiprocessing import Pool
from pathlib import Path
from itertools import product
from tqdm import tqdm



def _dock_one_wrapper(args):
    return _dock_one(*args)

def _dock_one(cfg, receptor, ligand, autobox_ligand, out_file, log_file):
    cmd = [
        cfg["paths"]["smina_path"],
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--autobox_ligand", str(autobox_ligand),
        "--autobox_add", str(cfg["docking_params"]["autobox_add"]),
        "--exhaustiveness", str(cfg["docking_params"]["exhaustiveness"]),
        "--num_modes", str(cfg["docking_params"]["num_modes"]),
        "--cpu", str(cfg["docking_params"]["cpu"]),
        "--out", str(out_file),
        "--log", str(log_file)
    ]

    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)



def _build_tasks(cfg):
    mode = cfg["docking_mode"]
    rec_dir = Path(cfg["paths"]["receptors_cleaned_folder"])
    lig_dir = Path(cfg["paths"]["ligands_folder"])
    nat_dir = Path(cfg["paths"]["native_ligands_folder"])
    out_dir = Path(cfg["paths"]["output_folder"])

    recs = sorted(rec_dir.glob("*.pdbqt"))
    ligs = sorted(lig_dir.glob("*.pdbqt"))

    nats_by_pdbid = {}
    for nat in nat_dir.glob("*.pdbqt"):
        parts = nat.stem.split("_")
        if len(parts) >= 2:
            pdbid = parts[0]
            nats_by_pdbid[pdbid] = nat

    tasks = []

    if mode == "redock_native":
        for pdbid, native in nats_by_pdbid.items():
            rec = rec_dir / f"{pdbid}.pdbqt"
            lig_id = native.stem.replace(f"{pdbid}_", "")
            if rec.exists():
                out = out_dir / f"{pdbid}__{lig_id}__native_redock.pdbqt"
                log = out_dir / f"{pdbid}__{lig_id}__native_redock.log"
                tasks.append((cfg, rec, native, native, out, log))

    elif mode == "diagonal":
        for lig in ligs:
            pdbid = lig.stem.split("_")[0]
            lig_id = lig.stem
            rec = rec_dir / f"{pdbid}.pdbqt"
            native = nats_by_pdbid.get(pdbid)
            if rec.exists() and native:
                out = out_dir / f"{pdbid}__{lig_id}__test_dock.pdbqt"
                log = out_dir / f"{pdbid}__{lig_id}__test_dock.log"
                tasks.append((cfg, rec, lig, native, out, log))

    elif mode == "full_matrix":
        for rec in recs:
            pdbid = rec.stem
            native = nats_by_pdbid.get(pdbid)
            if not native:
                continue
            for lig in ligs:
                lig_id = lig.stem
                out = out_dir / f"{pdbid}__{lig_id}__test_dock__.pdbqt"
                log = out_dir / f"{pdbid}__{lig_id}__test_dock__.log"
                tasks.append((cfg, rec, lig, native, out, log))

    return tasks



def run_batch_docking(cfg, log):
    tasks = _build_tasks(cfg)
    log.info("Dokowanie %d kompleksów...", len(tasks))

    for args in tqdm(tasks):
        try:
            print(f"[DOCKING] {args[1].name} vs {args[2].name} ...")  # receptor vs ligand

            _dock_one(*args)
        except Exception as e:
            print(f"❌ Error in docking {args}: {e}")


