from pathlib import Path
import pandas as pd
import pyranges as pr


gr_tes = pr.PyRanges(pd.read_csv(snakemake.input['tes']))
lapa_dir = Path(snakemake.input['lapa_dir'])

gr_neg = pr.read_bigwig(str(lapa_dir / 'all_tes_counts_counts_neg.bw')).df
gr_neg['Strand'] = '-'

gr_pos = pr.read_bigwig(str(lapa_dir / 'all_tes_counts_counts_pos.bw')).df
gr_pos['Strand'] = '+'

gr = pr.PyRanges(pd.concat([gr_pos, gr_neg]))

df = gr_tes.nearest(gr, strandedness='same').df
df['dist'] = df['End'] - df['End_b']
df = df[df['dist'].abs() < 50]

df_hist = df.groupby('dist')['Value'].sum().reset_index()

df_hist.reset_index().to_csv(snakemake.output['hist'], index=False)
with open(snakemake.output['percentage'], 'w') as f:
    f.write(str(df_hist['Value'].sum() / gr.Value.sum()))
