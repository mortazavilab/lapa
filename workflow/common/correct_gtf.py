from lapa.utils.io import read_talon_read_annot
from lapa.result import LapaResult
from lapa.read import read_tes_mapping, correct_gtf_tes, \
    tss_cluster, tss_mapping


print('PolyA site read mapping...')
df_cluster = LapaResult(snakemake.input['lapa_dir']).read_cluster()
df_mapping = read_tes_mapping(df_cluster,
                              snakemake.input['read_annot'])

print('TSS read mapping...')
df_tss_cluster = tss_cluster(snakemake.input['read_annot'],
                             snakemake.input['fasta'])
df_tss_mapping = tss_mapping(df_tss_cluster,
                             snakemake.input['read_annot'])

df_reads = read_talon_read_annot(snakemake.input['read_annot']).rename(
    columns={'annot_transcript_id': 'transcript_id'})

print('Correcting gtf...')
correct_gtf_tes(
    df_mapping,
    df_reads[['read_name', 'transcript_id']],
    snakemake.input['gtf'],
    snakemake.output['gtf'],
    df_read_tss=df_tss_mapping
)
