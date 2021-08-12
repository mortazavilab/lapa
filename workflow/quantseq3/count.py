from pathlib import Path
import pandas as pd


rows = list()

for path in snakemake.input['count_raw']:
    with open(path) as f:
        rows.append({
            'type': 'raw',
            'count': int(f.read())
        })

for path in snakemake.input['count_polyA']:
    with open(path) as f:
        rows.append({
            'type': 'polyA',
            'count': int(f.read())
        })

pd.DataFrame(rows).to_csv(snakemake.output['csv'], index=False)
