import pyranges as pr


def get_tes_from_gtf(gtf):
    gr_gtf = pr.read_gtf(gtf)

    df_tes = gr_gtf.features.tes().df[[
        'Chromosome', 'Start', 'End', 'Strand', 'gene_id', 'gene_type'
    ]].drop_duplicates()

    return df_tes
