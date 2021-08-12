from pathlib import Path
import pandas as pd


rows = list()

for path in snakemake.input['count_raw_dir']:
    with open(path) as f:
        rows.append({
            'encode_id': Path(path).stem,
            'type': 'raw',
            'count': int(f.read())
        })

for path in snakemake.input['count_dir']:
    with open(path) as f:
        rows.append({
            'encode_id': Path(path).stem,
            'type': 'polyA',
            'count': int(f.read())
        })

pd.DataFrame(rows).to_csv(snakemake.output['csv'], index=False)
