from lapa.utils.gtf import get_tes_from_gtf

df_tes = get_tes_from_gtf(snakemake.input['gtf'])
df_tes.to_csv(snakemake.output['tes'], index=False)
