from pathlib import Path
from itertools import product
import subprocess
from tqdm import tqdm

# ────────────────────────────────────────────────────────────────────
#  Param helpers
# ────────────────────────────────────────────────────────────────────
def _resolved_params(cfg) -> dict:
    """Zwraca słownik parametrów po połączeniu default + tryb."""
    mode = cfg["docking_mode"]
    params = dict(cfg["docking_params"].get("default", {}))   # kopia
    params.update(cfg["docking_params"].get(mode, {}))        # nadpisz
    return params


def _add_arg(cmd: list, flag: str, value):
    """Dodaj parę 'flag value' do cmd tylko jeśli value nie jest None."""
    if value is not None:
        cmd += [flag, str(value)]


# ────────────────────────────────────────────────────────────────────
#  Dokowanie pojedynczej pary
# ────────────────────────────────────────────────────────────────────
def _dock_one(cfg, receptor, ligand, autobox_ligand, out_file, log_file):
    params = _resolved_params(cfg)

    cmd = [
        cfg["paths"]["smina_path"],
        "--receptor",      str(receptor),
        "--ligand",        str(ligand),
        "--autobox_ligand", str(autobox_ligand),
        "--out",           str(out_file),
        "--log",           str(log_file),
    ]

    # opcjonalne argumenty – dorzucamy tylko jeśli są w YAML‑u
    _add_arg(cmd, "--autobox_add",  params.get("autobox_add"))
    _add_arg(cmd, "--exhaustiveness", params.get("exhaustiveness"))
    _add_arg(cmd, "--num_modes",      params.get("num_modes"))
    _add_arg(cmd, "--cpu",            params.get("cpu", cfg["runtime"].get("n_jobs", 1)))

    subprocess.run(cmd, check=True,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)


# ────────────────────────────────────────────────────────────────────
#  Budowanie listy zadań
# ────────────────────────────────────────────────────────────────────
def _build_tasks(cfg):
    mode     = cfg["docking_mode"]
    rec_dir  = Path(cfg["paths"]["receptors_cleaned_folder"])
    lig_dir  = Path(cfg["paths"]["active_ligands_folder"])
    nat_dir  = Path(cfg["paths"]["native_ligands_folder"])
    out_dir  = Path(cfg["paths"]["output_folder"])

    recs = sorted(rec_dir.glob("*.pdbqt"))
    ligs = sorted(lig_dir.glob("*.pdbqt"))

    # mapowanie PDB → ligand natywny
    native_map = {p.stem.split("_")[0]: p for p in nat_dir.glob("*.pdbqt")}

    tasks = []

    if mode == "redock_native":
        for pdb, nat in native_map.items():
            rec = rec_dir / f"{pdb}.pdbqt"
            if not rec.exists():
                continue
            lig_id = nat.stem.replace(f"{pdb}_", "")
            out = out_dir / f"{pdb}__{lig_id}__native_redock.pdbqt"
            log = out.with_suffix(".log")
            tasks.append((cfg, rec, nat, nat, out, log))

    elif mode == "diagonal":
        for lig in ligs:
            pdb = lig.stem.split("_")[0]
            rec = rec_dir / f"{pdb}.pdbqt"
            nat = native_map.get(pdb)
            if rec.exists() and nat:
                out = out_dir / f"{pdb}__{lig.stem}__diag.pdbqt"
                log = out.with_suffix(".log")
                tasks.append((cfg, rec, lig, nat, out, log))

    elif mode in {"matrix", "full_matrix"}:
        for rec, lig in product(recs, ligs):
            pdb = rec.stem
            nat = native_map.get(pdb)
            if not nat:
                continue
            out = out_dir / f"{pdb}__{lig.stem}__matrix.pdbqt"
            log = out.with_suffix(".log")
            tasks.append((cfg, rec, lig, nat, out, log))

    return tasks


# ────────────────────────────────────────────────────────────────────
#  Batch runner
# ────────────────────────────────────────────────────────────────────
def run_batch_docking(cfg, log):
    tasks = _build_tasks(cfg)
    log.info("Dokowanie %d kompleksów…", len(tasks))

    for args in tqdm(tasks):
        rec_name = args[1].name
        lig_name = args[2].name
        try:
            print(f"[DOCKING] {rec_name} vs {lig_name} …")
            _dock_one(*args)
        except Exception as e:
            print(f"❌  Error docking {rec_name} vs {lig_name}: {e}")
