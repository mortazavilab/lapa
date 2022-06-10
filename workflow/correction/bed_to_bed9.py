import pandas as pd


df = pd.read_csv(snakemake.input['bed'], sep='\t', header=None)
df[7] = 0
df[8] = 0
df.to_csv(snakemake.output['bed'], sep='\t', header=None, index=False)
