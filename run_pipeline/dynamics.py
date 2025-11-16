#!/usr/bin/env python3
"""
Minimal export of docking results for molecular dynamics.

This script builds a small, self‑contained directory with:
    - one native ligand per receptor (from `better_than_native.csv`)
    - the corresponding docked poses that outperformed the native

For each entry it stores the protein PDB, the ligand SDF and the SMILES
string in a `dynamics_manifest.csv` file, which can be consumed by
downstream MD workflows.
"""

import logging, shutil
from pathlib import Path
import pandas as pd
from openbabel import pybel

log = logging.getLogger(__name__)


# ───────────────── helpers ────────────────────────────────
def pdbqt_to_sdf(pdbqt: Path, sdf: Path) -> None:
    mol = next(pybel.readfile("pdbqt", str(pdbqt)))
    mol.addh(); mol.localopt(forcefield="uff", steps=500)
    mol.write("sdf", str(sdf), overwrite=True)


def normalize(df: pd.DataFrame) -> None:
    df.columns = [c.strip().lstrip("\ufeff") for c in df.columns]


# ───────────────── main ───────────────────────────────────
def prepare_for_dynamics(
    better_csv: Path,
    smiles_csv: Path,
    out_root: Path,
    visuals_root: Path,
    dynamics_root: Path,
    id_col="ID",
    smiles_col="SMILES",
) -> None:

    # 1) input tables -------------------------------------------------------
    better = pd.read_csv(better_csv)
    smiles = pd.read_csv(smiles_csv, encoding="utf-8-sig", sep=";", engine="python")
    normalize(smiles)

    id_col  = next(c for c in smiles if c.lower() == id_col.lower())
    smiles_col = next(c for c in smiles if c.lower() == smiles_col.lower())
    smiles_map = dict(zip(smiles[id_col], smiles[smiles_col]))

    # 2) output dirs (wipe ligands/ first) ----------------------------------
    base     = dynamics_root / "files_for_dynamics"
    lig_dir  = base / "ligands"
    prot_dir = base / "protein_cleaned"
    shutil.rmtree(lig_dir,  ignore_errors=True)
    shutil.rmtree(prot_dir, ignore_errors=True)
    lig_dir.mkdir(parents=True, exist_ok=True)
    prot_dir.mkdir(parents=True, exist_ok=True)
    rel = lambda p: p.relative_to(base)

    folder_map = {"native_redock": "dock_native",
                  "diag":          "diagonal",
                  "dock":          "matrix"}

    records = []

    # 3) one native ligand per receptor ------------------------------------
    for rec in better["receptor"].unique():
        nat = next((out_root / "dock_native").glob(f"{rec}__*__native_redock.pdbqt"), None)
        if not nat:
            continue
        lig_id = nat.stem.split("__")[1]
        sdf = lig_dir / f"{rec}__{lig_id}__native_redock.sdf"
        pdbqt_to_sdf(nat, sdf)

        vis = visuals_root / f"native_redock__{rec}__{lig_id}"
        rec_src = vis / f"{rec}_receptor.pdb"
        rec_dst = prot_dir / f"{rec}.pdb"
        if rec_src.exists():
            shutil.copy2(rec_src, rec_dst)

        records.append({
            "protein_pdb": str(rel(rec_dst)),
            "ligand_id":   lig_id,
            "ligand_sdf":  str(rel(sdf)),
            "smiles":      ""
        })

    # 4) better‑than‑native rows -------------------------------------------
    for _, row in better.iterrows():
        rec, lig, mode = row["receptor"], row["ligand"], row["mode"]
        sub = folder_map.get(mode, mode)
        pdbqt = out_root / sub / f"{rec}__{lig}__{mode}.pdbqt"
        if not pdbqt.exists():
            continue
        sdf = lig_dir / f"{rec}__{lig}__{mode}.sdf"
        pdbqt_to_sdf(pdbqt, sdf)

        vis  = visuals_root / f"{mode}__{rec}__{lig}"
        rec_src = vis / f"{rec}_receptor.pdb"
        rec_dst = prot_dir / f"{rec}.pdb"
        if rec_src.exists() and not rec_dst.exists():
            shutil.copy2(rec_src, rec_dst)

        records.append({
            "protein_pdb": str(rel(rec_dst)),
            "ligand_id":   lig,
            "ligand_sdf":  str(rel(sdf)),
            "smiles":      smiles_map.get(lig, "")
        })

    # 5) manifest -----------------------------------------------------------
    pd.DataFrame(records,
                 columns=["protein_pdb", "ligand_id", "ligand_sdf", "smiles"]
                 ).to_csv(base / "dynamics_manifest.csv", index=False)
    log.info("Manifest ready (%d entries)", len(records))


def _default_paths() -> tuple[Path, Path, Path, Path, Path]:
    """
    Construct default locations relative to the repository root.

    The layout assumes this file lives in `run_pipeline/` and that the
    standard folders (`data/`, `results/`) follow the structure described
    in the project README.
    """
    root = Path(__file__).resolve().parent
    better_csv = root / "results" / "better_than_native.csv"
    smiles_csv = root / "data" / "ligands" / "ligands.csv"
    out_root = root / "results"
    visuals_root = root / "results" / "visuals"
    dynamics_root = root
    return better_csv, smiles_csv, out_root, visuals_root, dynamics_root


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    better_csv, smiles_csv, out_root, visuals_root, dynamics_root = _default_paths()

    prepare_for_dynamics(
        better_csv,
        smiles_csv,
        out_root,
        visuals_root,
        dynamics_root,
        id_col="ID",
        smiles_col="SMILES",
    )
