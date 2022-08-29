import time
import pandas as pd

df = list()

for cell, encode_ids in snakemake.config['c2c12']['bam'].items():
    for encode_id in encode_ids:
        df.append({
            'sample': encode_id,
            'dataset': cell,
            'path': snakemake.config['encode']['bam'].format(encode_id=encode_id)
        })

pd.DataFrame(df).to_csv(snakemake.output['sample_config'], index=False)

time.sleep(3) # needed for snakemake
