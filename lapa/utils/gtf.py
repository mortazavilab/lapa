import pyranges as pr


def get_tes_from_gtf(gtf):
    '''
    Extract unique transcript ends sites from gtf files as datafarme.

    Args:
      gtf: Path to gtf file.
    '''
    gr_gtf = pr.read_gtf(gtf)

    df_tes = gr_gtf.features.tes().df[[
        'Chromosome', 'Start', 'End', 'Strand', 'gene_id', 'gene_type'
    ]].drop_duplicates()

    return df_tes


def get_tss_from_gtf(gtf):
    '''
    Extract unique transcript start sites from gtf files as datafarme.

    Args:
      gtf: Path to gtf file.
    '''
    gr_gtf = pr.read_gtf(gtf)

    df_tss = gr_gtf.features.tss().df[[
        'Chromosome', 'Start', 'End', 'Strand', 'gene_id', 'gene_type'
    ]].drop_duplicates()

    return df_tss
