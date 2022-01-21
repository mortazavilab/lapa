import pdb
import pyranges as pr


df_gtf = pr.read_gtf(snakemake.input['gtf']).df
df_gtf = df_gtf[df_gtf['Feature'] == 'gene']

df_gtf[['gene_id', 'gene_name']].to_csv(
    snakemake.output['mapping'], index=False)
