# Molecular Docking Pipeline

> \*\*Fast, idempotent & fully‑scriptable molecular docking workflow built around [Smina](https://github.com/mwojcikowski/smina) + AutoDockTools (py3).
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
  git clone https://github.com/<you>/molecular-docking.git
  cd molecular-docking

  # 1. Conda is recommended
  conda env create -f environment.yaml
  conda activate docking

  # 2. Over‑the‑air tools
  # (edit paths in config.yaml if you install elsewhere)
  pip install smina
  # GUI PyMOL (on macOS)
  brew install --cask pymol
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
├── config.yaml
├── data/
│   ├── crystals/           # raw .pdb from RCSB
│   ├── receptors/          # cleaned PDB
│   ├── receptors_qt/       # PDBQT ready for Smina
│   ├── ligands/            # test ligands (PDBQT) computed from smiles csv
│   └── native_ligands/     # native co‑crystallised ligands (PDBQT)
├── results/
│   ├── dock_native/        # *.pdbqt + .log for each dock
│   ├── matrix/
│   ├── visuals/            # RMSD plots, SVG overlays, Py3Dmol files
│   ├── results.csv         # consolidated docking table 
│   └── 3D_visuals          # Jupiter notebook for interactive docking results visualizations 
└── run.log
```

---

## ⚙️ Configuration (`config.yaml`)

```yaml
docking_mode: redock_native   # or diagonal/matrix/full_matrix

ligand_csv:
  enabled: true               # set true when you have a CSV
  file: data/ligands/smiles.csv
  smiles_column: SMILES
  id_column: ID

paths:
  receptors_folder: data/receptors
  receptors_cleaned_folder: data/receptors_qt
  crystals_folder: data/crystals
  ligands_folder: data/ligands
  native_ligands_folder: data/natives
  output_folder: results
  smina_path: /usr/local/bin/smina
  adt_root:  /opt/AutoDockTools
  pymol_path: /Applications/PyMOL.app/Contents/MacOS/PyMOL
  
docking_params:
  default:
    exhaustiveness: 8
    num_modes: 9
  redock_native:
    autobox_add: 4
    
runtime:
  n_jobs: 4                  # CPU threads for Smina

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
| **`utils.file_utils`**                | `ensure_dirs`, `setup_logging`         | Small helpers for filesystem & logging.                            |
| **`utils.validation`**                | `files_ok`, `rmsd_atoms`               | Lightweight success checks.                                        |
| **`run_pipeline.main`**               | `main()`                               | Glue everything together.                                          |

---




## 🤝 Contributing

Pull requests are welcome! Please open an issue first to discuss major changes.
Make sure to run `black .` and `ruff --fix .` before committing.

---

## 📣 Citation

If you find this pipeline useful in your research, please cite:

```
@software{docking_pipeline_2025,
  author  = Tomasz Szostek,
  title   = Molecular docking pipeline,
  year    = 2025,
  url     = {https://github.com/TomaszSzostek/Molecular_docking.git}
}
```

---

*Happy docking! 🔬✨*
