import numpy as np
import pandas as pd
import pyranges as pr
from lapa.result import LapaResult


gr_tes = pr.read_gtf(snakemake.input['gtf']).features.tes()

result = LapaResult(snakemake.input['lapa_dir'])

df_quantseq = result.read_clusters().drop_duplicates(
    subset=['Chromosome', 'polyA_site', 'Strand'])
df_quantseq['Start'] = df_quantseq['polyA_site'] - 1
df_quantseq['End'] = df_quantseq['polyA_site']

# Find overlap between TES and quantseq polyA cluster
df_overlap = gr_tes.nearest(pr.PyRanges(df_quantseq)).df
dist = snakemake.params['max_distance']
df_overlap = df_overlap[df_overlap['Distance'] < dist]

# df_overlap = df_overlap[df_overlap['TPM_average'] >
#                         float(snakemake.wildcards['tpm_filter'])]

transcript_overlap = set(df_overlap['transcript_id'])

# Prepare output file which contains transcript_id, support, novelty
df_abundance = pd.read_csv(snakemake.input['abundance'], sep='\t') \
                 .rename(columns={'annot_transcript_id': '_transcript_id'}) \
                 .set_index('_transcript_id')

df = gr_tes.df
df['_transcript_id'] = df['transcript_id'].str.split('#').str.get(0)
df = df.set_index('_transcript_id').join(df_abundance, how='inner')

df.loc[df['ISM_subtype'].isin(
    {'None', 'Both'}), 'ISM_subtype'] = 'other'
df['ISM_subtype'] = df['ISM_subtype'].str.lower()

df['transcript_novelty'] = np.where(
    df['transcript_novelty'] == 'ISM',
    df['ISM_subtype'] + ' ' + df['transcript_novelty'],
    df['transcript_novelty']
)

df = df.set_index('transcript_id')
transcripts = df.index

# overlap plot
pd.DataFrame({
    'transcript_id': transcripts,
    'support': np.where(transcripts.isin(transcript_overlap), 'yes', 'no'),
    'novelty': df.loc[transcripts, 'transcript_novelty'].tolist()
}).to_csv(snakemake.output['overlap'], index=False)
