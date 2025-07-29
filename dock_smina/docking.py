"""
dock_smina/docking.py

Build and execute batch docking tasks with Smina according to the configured mode.

This module provides:
  - Resolving and merging default and mode‑specific docking parameters.
  - Single‑ligand docking invocation with proper Smina command‑line flags.
  - Generation of task lists for 'redock_native', 'diagonal' and 'full_matrix' modes.
  - Sequential execution of all tasks with a progress bar.

Public API:
    run_batch_docking(cfg: dict, log: logging.Logger) -> None

Configuration keys used:
    cfg["docking_mode"]
    cfg["docking_params"]
    cfg["runtime"]["n_jobs"]
    cfg["paths"]["smina_path"]
    cfg["paths"]["receptors_cleaned_folder"]
    cfg["paths"]["active_ligands_folder"]
    cfg["paths"]["native_ligands_folder"]
    cfg["paths"]["output_folder"]
"""

from pathlib import Path
from itertools import product
import subprocess
from tqdm import tqdm

def _resolved_params(cfg) -> dict:
    """
    Merge default and mode‑specific docking parameters.

    Parameters
    ----------
    cfg : dict
        Pipeline configuration containing:
          - 'docking_params.default': base parameters (dict)
          - 'docking_params[mode]': overrides for the current mode (dict)
          - 'docking_mode': name of the current mode (str)

    Returns
    -------
    dict
        Combined parameter set for Smina flags.
    """
    mode = cfg["docking_mode"]
    params = dict(cfg["docking_params"].get("default", {}))
    params.update(cfg["docking_params"].get(mode, {}))
    return params

def _add_arg(cmd: list, flag: str, value):
    """
    Append a flag and its value to the command list if the value is not None.

    Parameters
    ----------
    cmd : list
        List of command‑line arguments being constructed.
    flag : str
        Command‑line flag (e.g., '--exhaustiveness').
    value : any
        Value for the flag; if None, the flag is omitted.
    """
    if value is not None:
        cmd += [flag, str(value)]

def _dock_one(cfg, receptor, ligand, autobox_ligand, out_file, log_file):
    """
    Run a single Smina docking job with specified parameters.

    Parameters
    ----------
    cfg : dict
        Pipeline configuration (paths, runtime, docking_params).
    receptor : Path
        Path to the receptor PDBQT file.
    ligand : Path
        Path to the ligand PDBQT file.
    autobox_ligand : Path or None
        Reference ligand file for defining the docking box (if applicable).
    out_file : Path
        Output PDBQT path for the resulting docked pose.
    log_file : Path
        File path for Smina log output.

    Raises
    ------
    subprocess.CalledProcessError
        If Smina execution fails (due to check=True).
    """
    params = _resolved_params(cfg)

    cmd = [
        cfg["paths"]["smina_path"],
        "--receptor", str(receptor),
        "--ligand",   str(ligand),
        "--out",      str(out_file),
        "--log",      str(log_file),
    ]

    # only add the flag when we actually have a reference ligand
    if autobox_ligand is not None:
        cmd += ["--autobox_ligand", str(autobox_ligand)]

    # add optional docking parameters
    _add_arg(cmd, "--autobox_add",    params.get("autobox_add"))
    _add_arg(cmd, "--exhaustiveness", params.get("exhaustiveness"))
    _add_arg(cmd, "--num_modes",      params.get("num_modes"))
    _add_arg(cmd, "--cpu",            params.get("cpu", cfg["runtime"].get("n_jobs", 1)))

    # execute Smina silently
    subprocess.run(cmd, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def _build_tasks(cfg):
    """
    Construct the list of docking tasks based on the selected mode.

    Parameters
    ----------
    cfg : dict
        Pipeline configuration with keys:
          - 'docking_mode': current docking mode (str)
          - 'paths.receptors_cleaned_folder': prepared receptor folder
          - 'paths.active_ligands_folder': test ligands folder
          - 'paths.native_ligands_folder': native ligands folder
          - 'paths.output_folder': base output directory

    Returns
    -------
    list of tuple
        Each tuple contains arguments for a single docking job:
        (cfg, receptor, ligand, autobox_ligand, out_file, log_file).
    """
    mode     = cfg["docking_mode"]
    rec_dir  = Path(cfg["paths"]["receptors_cleaned_folder"])
    lig_dir  = Path(cfg["paths"]["active_ligands_folder"])
    nat_dir  = Path(cfg["paths"]["native_ligands_folder"])
    out_dir  = Path(cfg["paths"]["output_folder"])

    recs = sorted(rec_dir.glob("*.pdbqt"))
    ligs = sorted(lig_dir.glob("*.pdbqt"))

    # map receptor ID to native ligand file (based on filename prefix)
    native_map = {p.stem.split("_")[0]: p for p in nat_dir.glob("*.pdbqt")}
    tasks = []

    if mode == "redock_native":
        subdir = "dock_native"
        task_dir = out_dir / subdir
        task_dir.mkdir(parents=True, exist_ok=True)
        for pdb, nat in native_map.items():
            rec = rec_dir / f"{pdb}.pdbqt"
            if not rec.exists():
                continue
            lig_id = nat.stem.replace(f"{pdb}_", "")
            out = task_dir / f"{pdb}__{lig_id}__native_redock.pdbqt"
            log = out.with_suffix(".log")
            if out.exists():
                continue
            tasks.append((cfg, rec, nat, nat, out, log))

    elif mode == "diagonal":
        subdir = "diagonal"
        task_dir = out_dir / subdir
        task_dir.mkdir(parents=True, exist_ok=True)
        for lig in ligs:
            pdb = lig.stem.split("_")[0]
            rec = rec_dir / f"{pdb}.pdbqt"
            nat = native_map.get(pdb)
            if rec.exists() and nat:
                out = task_dir / f"{pdb}__{lig.stem}__diag.pdbqt"
                log = out.with_suffix(".log")
                if out.exists():
                    continue
                tasks.append((cfg, rec, lig, nat, out, log))

    elif mode in {"matrix", "full_matrix"}:
        subdir = "matrix"
        task_dir = out_dir / subdir
        task_dir.mkdir(parents=True, exist_ok=True)
        for rec, lig in product(recs, ligs):
            pdb = rec.stem
            nat = native_map.get(pdb)
            out = task_dir / f"{pdb}__{lig.stem}__dock.pdbqt"
            log = out.with_suffix(".log")
            if out.exists():
                continue
            tasks.append((cfg, rec, lig, nat, out, log))

    return tasks

def run_batch_docking(cfg, log):
    """
    Discover and execute all docking tasks for the current mode.

    Parameters
    ----------
    cfg : dict
        Pipeline configuration.
    log : logging.Logger
        Logger for informational and error messages.

    Returns
    -------
    None
    """
    log.info(">>> DOCK")
    tasks = _build_tasks(cfg)
    if not tasks:
        log.info("Docking already performed– skipping...")
        return

    # iterate through tasks with progress bar
    for task in tqdm(tasks, desc="Docking"):
        try:
            _dock_one(*task)
        except Exception as e:
            log.error(f"Error during docking: {e}")

