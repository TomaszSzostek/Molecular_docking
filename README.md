# Molecular Docking Pipeline

> **Fast**, idempotent & fully‑scriptable molecular docking workflow built around [Smina](https://github.com/mwojcikowski/smina) + AutoDockTools (py3).
> From raw PDB IDs & SMILES to ranked hits, RMSD overlays and publication‑ready interactive visualisations — in **one command**.

---

## ✨ Highlights

| Feature                   | What it does                                                                                                                            |
| ------------------------- |-----------------------------------------------------------------------------------------------------------------------------------------|
| **End‑to‑end automation** | Fetches crystal structures, prepares PDBQT, runs batch docking (3 modes), consolidates logs, calculates RMSD and auto‑generates visuals |
| **Idempotent by design**  | Re‑running the pipeline never repeats an already‑finished step — handy for large virtual screens                                        |
| **Config‑driven**         | Everything lives in a single `config.yaml`: paths, run modes, docking params, CPU count…                                                |
| **Lean dependencies**     | Pure‑Python except for Smina executables;                                                                                               |
| **Reproducible results**  | Dockerfile & Conda env supplied; every output artefact is timestamped and logged                                                        |

---

## 📦 Quick install

```bash
git clone https://github.com/TomaszSzostek/Molecular_docking.git
cd Molecular_docking

# 1. Create Conda environment
conda env create -f environment.yml
conda activate docking

# 2. Install external tools (optional, system‑wide)
# Smina binaries (or use your own build)
#   https://github.com/mwojcikowski/smina
# PyMOL for 3D inspection (example for macOS)
#   brew install --cask pymol
```


---

## 🚀 One‑liner usage

```bash
  python run_pipeline/main.py --config config.yaml
  # add --dry-run to only print planned actions
```

By default the pipeline runs in **redock\_native** mode (native ligand vs receptor).
Switch to **diagonal** or **full\_matrix** in the YAML to explore cross‑docking scenarios.

---

## 🗂️ Directory layout (after first run)

```
├── run_pipeline/config.yaml
├── run_pipeline/data/
│   ├── crystals/           # raw .pdb from RCSB
│   ├── receptors/          # cleaned PDB (PDB)
│   ├── receptors_cleaned/  # PDBQT ready for Smina
│   ├── ligands/            # test ligands (PDBQT) from ligands.csv
│   └── native_ligands/     # native co‑crystallised ligands (PDB/PDBQT)
├── run_pipeline/results/
│   ├── dock_native/        # *.pdbqt + .log for native redocks
│   ├── matrix/             # *.pdbqt + .log for full_matrix / diagonal runs
│   ├── visuals/            # RMSD plots, overlays, 2D boards
│   ├── results.csv         # consolidated docking table
│   ├── rmsd_summary.csv    # per‑receptor RMSD statistics
│   └── better_than_native.csv
└── run_pipeline/run.log
```

---

## ⚙️ Configuration (`run_pipeline/config.yaml`)

```yaml
docking_mode: redock_native   # or diagonal/full_matrix

ligand_csv:
  enabled: true               # set true when you have a CSV with test ligands (full_matrix/diagonal)
  file: "data/ligands/ligands.csv"
  smiles_column: "SMILES"
  id_column: "ID"

paths:
  receptors_folder: "data/receptors"
  receptors_cleaned_folder: "data/receptors_cleaned"
  native_ligands_folder: "data/native_ligands"
  ligands_folder: "data/ligands"
  output_folder: "results"
  crystals_folder: "data/crystals"
  temp_dir: "tmp"
  smina_path: "/usr/local/bin/smina"
  visuals: "results/visuals"

docking_params:
  default:
    num_modes: 6
    exhaustiveness: 16
    autobox_add: 1

  redock_native:
    num_modes: 12
    exhaustiveness: 16
    autobox_add: 4
    
runtime:
  pool_size: 6               # concurrent docking workers
  cpu_per_job: 1             # threads passed to Smina via --cpu
```

---

## 🛠️ Module‑by‑module tour

| Module                                | Public API                             | Purpose                                                            |
| ------------------------------------- | -------------------------------------- |--------------------------------------------------------------------|
| **`prepare_inputs.fetch_crystals`**   | `fetch_and_split_batch(cfg, log)`      | Download PDB structures and split receptor/ligand chains.          |
| **`prepare_inputs.adt_preparator`**   | `receptor_to_pdbqt`, `ligand_to_pdbqt` | Run AutoDockTools scripts + RDKit sanitisation to create PDBQT.    |
| **`prepare_inputs.smiles_loader`**    | `csv_to_pdbqt`                         | Convert a CSV (ID,SMILES) into 3‑D ligand PDBQT batch.             |
| **`dock_smina.docking`**              | `run_batch_docking`                    | Build docking task list and invoke Smina with correct params.      |
| **`analyze.results_extractor`**       | `consolidate_logs`                     | Parse all \*.log files → `results.csv` (min/avg binding energies). |
| **`analyze.postprocess`**             | `rank_vs_native`                       | Flag ligands that beat the native by an affinity margin.           |
| **`analyze.RMSD`**                    | `run_rmsd_and_plot`                    | Kabsch RMSD vs native + 2‑D overlay (RDKit).                       |
| **`analyze.files_for_visualization`** | `generate_files`                       | Prep Py3Dmol & Jupyter‑friendly files for manual inspection.       |
| **`visualize_2d`**                    | `python -m visualize_2d --complex …`   | Pastel 2D PoseView‑style protein–ligand interaction boards.        |
| **`utils.file_utils`**                | `ensure_dirs`, `setup_logging`         | Small helpers for filesystem & logging.                            |
| **`utils.validation`**                | `files_ok`, `rmsd_atoms`               | Lightweight success checks.                                        |
| **`run_pipeline.main`**               | `main()`                               | Glue everything together.                                          |

### ✨ Pastel 2D interaction boards

After running the pipeline (which produces `results/visuals/complex__<mode>__...` folders), you can render custom 2D diagrams:

```
python -m visualize_2d \
  --complex run_pipeline/results/visuals/complex__matrix__1M17__1 \
  --out run_pipeline/results/visuals/2d_boards/1M17__1.png
```

The renderer blends RDKit stick figures, PLIP interactions and a pastel UI kit. Outputs can be `png` or `svg` and embed the auto-generated 3D snapshot if present.

---





## 📣 Citation

If you find this pipeline useful in your research, please cite:

```

  author  = Tomasz Szostek,
  title   = Molecular docking pipeline,
  year    = 2025,
  url     = https://github.com/TomaszSzostek/Molecular_docking.git

```

---

*Happy docking! 🔬✨*
