# 🔬 Docking Pipeline
A production‑ready, parallel **molecular‑docking workflow** (Smina + AutoDockTools) controlled by a single YAML file. Works on macOS (Apple Silicon/Intel) and Linux.

---
## 🔑 Key features
* **Single‑file configuration** – edit `config.yaml`, no code changes.
* **Pluggable AutoDockTools path** – set env var `ADT_HOME` or `paths.adt_home`.
* **Three docking modes** – `full_matrix`, `diagonal`, `redock_native`.
* **Automatic COM grid boxes** around native ligands.
* **Parallel execution** with live progress bar and auto‑retry of failed jobs.
* **Full results export** (all poses) + Top‑N summary 
* **Poses visualization and RMSD calculation** for redocked native ligands. 

---
## 🚀 Start
### 1 Prerequisites
* Linux/macOS/WSL (Windows only if Smina/Vina compiled).  
* Conda ≥ 4.10.
* Install Smina separately (e.g. Homebrew or system binary) and set smina_path in config.yaml”
* Install AutoDockTools separately (from available git repository)

## 📦 Installation
### 1 Clone AutoDockTools_py3 locally

```bash
  # choose any folder you like
git clone https://github.com/<your‑fork>/AutoDockTools_py3.git ~/tools/AutoDockTools_py3
export ADT_HOME=~/tools/AutoDockTools_py3   # add to ~/.zshrc for permanence
````
###Install Smina Separately 
```bash
  # macOS
brew install smina
  # Linux (x86_64)
conda install -c bioconda vina-smina
which smina   # note the path
```

### Create Conda environment 
```bash 
  conda env create -f environment.yml -n docking
  conda activate docking
  python - <<'PY'
  import rdkit, os; print('✓ RDKit OK, ADT_HOME=', os.getenv('ADT_HOME'))
PY
```
## ✨ Project Layout
```
project_root/
│
├── docking_pipeline/          
│   ├── __init__.py
│   ├── main.py
│   ├── config.yaml  
│   │── environment.yml
│   └── data/                    
│      ├── receptors/            
│      ├── native_ligands/        
│      ├── ligands/              
│      └── output/               
│          ├── *.pdbqt
│          ├── *.log
│          ├── results.parquet
│          ├── poses_min_max.csv
│          ├── hits_vs_native.csv
│          └── run.log
│   
├── utils/
│   ├── __init__.py
│   ├── file_utils.py
│   └── validation.py
│   
├── prepare/
│   ├── __init__.py
│   ├── adt_preparator.py  
│   ├── smlies_loader.py
│   └── fetch_crystals.py
│   
├── dock/
│   ├── __init__.py
│   └── docking.py
│   
├── analyze/
│   ├── __init__.py
│   ├── results_extractor.py
│   ├── RMSD.py
│   └── postprocess.py
│          
├── README.md

```

## 🧬 Quick run  
```bash
  cd project_root
  python -m docking_pipeline.main --config config.yaml
```
