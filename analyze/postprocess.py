import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def rank_vs_native(cfg, log):
    outdir = Path(cfg['paths']['output_folder'])
    df = pd.read_parquet(outdir / 'results.parquet')

    # native & hit tables
    native = df[df.is_native].groupby('receptor')['score'].min().rename('native_score')
    hits   = df[~df.is_native].merge(native, on='receptor')

    # save min/max per ligand × receptor  (full supplement)
    per_lig = df.groupby(['ligand','receptor']) \
                .agg(min_affinity=('score','min'),
                     max_affinity=('score','max')) \
                .reset_index()
    per_lig.to_csv(outdir / 'poses_min_max.csv', index=False)

    # filter by native_margin
    margin = cfg['postprocess']['native_margin']
    hits = hits[hits.score < hits.native_score - margin]

    # global Top‑N for the manuscript
    top_n = cfg['postprocess']['top_n']
    hits.nsmallest(top_n, 'score') \
        .to_csv(outdir / 'hits_vs_native.csv', index=False)

    # histogram
    plt.figure()
    hits['score'].hist(bins=30)
    plt.xlabel('Docking score (kcal/mol)')
    plt.ylabel('Count')
    plt.title('Scores better than native')
    plt.savefig(outdir / 'score_hist.png', dpi=300)
    plt.close()               # <‑‑ free memory

    log.info(f"Saved: poses_min_max.csv ({len(per_lig)} rows), "
             f"hits_vs_native.csv ({len(hits)}), score_hist.png")