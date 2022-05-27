import pandas as pd
from lapa.count import PolyaTailCounter


dists = [
    PolyaTailCounter(i, min_tail_len=1, count_aligned=True).tail_len_dist()
    for i in snakemake.input['bams']
]

df = pd.concat(
    pd.DataFrame({
        'len': dist.index.tolist(),
        'count': dist.tolist(),
    })
    for i, dist in enumerate(dists)
).groupby(['len']).sum().reset_index()

df['cumulative count'] = df['count'].cumsum()
df = df[df['cumulative count'] / df['count'].sum() < 0.99]

df['data_source'] = f'{snakemake.wildcards["platform"]} {snakemake.wildcards["library_prep"]}'

df.to_csv(snakemake.output['counts'], index=False)
