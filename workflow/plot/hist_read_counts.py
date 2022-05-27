import pandas as pd
import pyranges as pr
from pathlib import Path
from lapa.result import LapaResult


gr_tes = pr.PyRanges(pd.read_csv(snakemake.input['tes']))
gr = pr.PyRanges(LapaResult(snakemake.input['lapa_dir']).read_counts())

df = gr_tes.nearest(gr, strandedness='same').df
df['dist'] = df['End'] - df['End_b']
df = df[df['dist'].abs() < 50]

df_hist = df.groupby('dist')['count'].sum().reset_index()

Path(snakemake.output['hist']).parent \
                              .mkdir(exist_ok=True, parents=True)
df_hist.reset_index().to_csv(snakemake.output['hist'], index=False)

Path(snakemake.output['percentage']).parent \
                                    .mkdir(exist_ok=True, parents=True)
with open(snakemake.output['percentage'], 'w') as f:
    f.write(str(df_hist['count'].sum() / gr.count.sum()))
