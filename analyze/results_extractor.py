import re, pandas as pd
from pathlib import Path

aff = re.compile(r"REMARK VINA RESULT:\s+([\-\d\.]+)")


def parse(log: Path):
    txt = log.read_text()
    m = aff.search(txt)
    score = float(m.group(1)) if m else None
    parts = log.stem.split('__to__')
    lig, rec, tag = parts[0], parts[1], parts[2]
    return {'ligand':lig,'receptor':rec,'score':score,'is_native': tag=='native'}


def consolidate_logs(cfg, log):
    out = Path(cfg['paths']['output_folder'])
    logs = list(out.glob('*.log'))
    log.info(f'Parsing {len(logs)} logs …')
    recs = [parse(lg) for lg in logs]
    df = pd.DataFrame(recs)
    df.to_parquet(out/'results.parquet', index=False)
    df.to_csv(out/'results.csv',index=False)