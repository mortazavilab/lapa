import pyranges as pr

df_gtf = pr.read_gtf(snakemake.input['gtf'])

df_tes = df_gtf.features.tes().df[[
    'Chromosome', 'Start', 'End', 'Strand', 'gene_id', 'gene_type'
]].drop_duplicates().to_csv(snakemake.output['tes'], index=False)
