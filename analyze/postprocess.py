import sys, pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import SDWriter
from utils.validation import rmsd_atoms

def rank_vs_native(cfg, log):
    out = Path(cfg["paths"]["output_folder"])
    df_path = out / "results.csv"
    if not df_path.exists():
        log.error("results.csv not found. Run docking first.")
        return

    df = pd.read_csv(df_path)
    margin = cfg["postprocess"]["native_margin"]

    # tylko native redock
    if "redock" not in df["mode"].iloc[0]:
        log.info("Skipping postprocess – not native redocking mode.")
        return

    # znajdź minimalne affinitiy dla natywnego liganda (per receptor)
    native_df = df[df["is_native"] == True]
    native_scores = native_df.groupby("receptor")["min_affinity"].min().rename("native_score")

    # znajdź wszystkie testowe ligandy (czyli nie-native)
    hits = df[df["is_native"] == False].copy()

    # dołącz natywne affinitiy do każdej grupy receptorowej
    hits = hits.merge(native_scores, on="receptor", how="left")

    # filtruj tylko te, które mają lepsze wiązanie niż natywny ligand - z marginesem
    hits = hits[hits["min_affinity"] < hits["native_score"] - margin]

    hits_path = out / "better_than_native.csv"
    hits.to_csv(hits_path, index=False)
    log.info(f"Found {len(hits)} hits better than native ligands.")
