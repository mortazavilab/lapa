import pandas as pd
from lapa.result import LapaResult


result = LapaResult(snakemake.input['long_read'])

df = result.fisher_exact_test({
    'undif': ['ENCFF772LYG', 'ENCFF421MIL'],
    'dif': ['ENCFF699KOR', 'ENCFF731HHB']
}, min_gene_count=50)

df.reset_index() \
  .rename(columns={'index': 'polya_site'}) \
  .merge(pd.read_csv(snakemake.input['mapping']), on='gene_id') \
  .to_csv(snakemake.output['stats'], index=False)
