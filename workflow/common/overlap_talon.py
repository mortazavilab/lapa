import pdb
import numpy as np
import pandas as pd
import pyranges as pr
from lapa.utils.gtf import get_tes_from_gtf
from lapa.result import LapaResult


gr_tes = pr.read_gtf(snakemake.input['gtf']).features.tes()

result = LapaResult(snakemake.input['lapa_dir'])

df_quantseq = pd.concat([
    result.read_apa(i)
    for i in snakemake.params['samples']
]).drop_duplicates(subset=['Chromosome', 'polyA_site', 'Strand'])

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
                 .rename(columns={'annot_transcript_id': 'transcript_id'}) \
                 .set_index('transcript_id')

df_abundance.loc[df_abundance['ISM_subtype'].isin(
    {'None', 'Both'}), 'ISM_subtype'] = 'other'
df_abundance['ISM_subtype'] = df_abundance['ISM_subtype'].str.lower()

df_abundance['transcript_novelty'] = np.where(
    df_abundance['transcript_novelty'] == 'ISM',
    df_abundance['ISM_subtype'] + ' ' + df_abundance['transcript_novelty'],
    df_abundance['transcript_novelty']
)

transcripts = df_abundance.index.intersection(gr_tes.transcript_id)

# overlap plot
pd.DataFrame({
    'transcript_id': transcripts,
    'support': np.where(transcripts.isin(
        transcript_overlap), 'yes', 'no'),
    'novelty': df_abundance.loc[
        transcripts, 'transcript_novelty'].tolist()
}).to_csv(snakemake.output['overlap'], index=False)
