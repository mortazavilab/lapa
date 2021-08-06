import pandas as pd

df = pd.read_csv(snakemake.input['read_annot'], sep='\t')

bulk = ['PB155', 'PB213', 'PB214', 'PB154']

df[df['dataset'].isin(bulk)].to_csv(
    snakemake.output['bulk_annot'], sep='\t', index=False)
df[~df['dataset'].isin(bulk)].to_csv(
    snakemake.output['sc_annot'], sep='\t', index=False)
