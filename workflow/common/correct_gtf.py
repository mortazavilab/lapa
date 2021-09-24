import pdb
import pandas as pd
from lapa.utils.io import read_talon_read_annot
from lapa.read import read_tes_mapping, correct_gtf_tes


print('Reading input...')
df_mapping = read_tes_mapping(
    snakemake.input['lapa_dir'],
    snakemake.input['read_annot']
)
df_reads = read_talon_read_annot(snakemake.input['read_annot']).rename(
    columns={'annot_transcript_id': 'transcript_id'})

# pd.read_csv(snakemake.input['abundance'], sep='\t').rename(
#    columns={'annot_transcript_id': 'transcript_id'})

print('Correcting gtf...')
correct_gtf_tes(
    df_mapping,
    df_reads[['read_name', 'transcript_id']],
    snakemake.input['gtf'],
    snakemake.output['gtf']
)
