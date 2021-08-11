import numpy as np
import pandas as pd
import pyranges as pr


_core_columns = ['Chromosome', 'Start', 'End', 'Strand']


class GenomicRegions:

    def __init__(self, gtf_file):
        if type(gtf_file) == str:
            self.gr = pr.read_gtf(gtf_file)
        elif type(gtf_file) == pr.PyRanges:
            self.gr = gtf_file
        else:
            raise ValueError('invalid gtf_file argument')

    def features(self, features: set = None) -> pr.PyRanges:
        features = features or {'five_prime_utr', 'exon',
                                'gene', 'three_prime_utr', 'UTR'}
        gr = self.gr[self.gr.Feature.isin(features)]
        gr_introns = self.gr.features.introns()

        df = pd.concat([gr.df, gr_introns.df])
        df = df[[*_core_columns, 'Feature', 'gene_id', 'gene_name']]

        df['canonical_site'] = np.where(
            df['Feature'] == 'three_prime_utr',
            np.where(df['Strand'] == '-', df['Start'], df['End']),
            -1)

        return df

    def annotate(self, gr, features=None):
        df_gtf = self.features(features)
        df_gtf = df_gtf[df_gtf['Strand'].isin(gr.Strand.cat.categories)]
        gr_gtf = pr.PyRanges(df_gtf, int64=True)

        gr_ann = pr.PyRanges(gr.df, int64=True).join(
            gr_gtf, strandedness='same', how='left') \
            .drop(['Start_b', 'End_b', 'Strand_b'])
        df = gr_ann.df

        # need to avoid duplication of groupby for category
        df['Chromosome'] = df['Chromosome'].astype(str)
        df['Strand'] = df['Strand'].astype(str)

        df['Feature'] = df['Feature'].replace('-1', 'intergenic')
        df['gene_id'] = df['gene_id'].replace('-1', '')
        df['gene_name'] = df['gene_name'].replace('-1', '')

        # if gene_id defined overlap
        df = df.groupby(_core_columns).apply(self._agg_annotation_gene) \
                                      .reset_index(drop=True)
        return df

    @staticmethod
    def _agg_annotation_gene(df):
        feature_order = ['three_prime_utr',  'UTR', 'exon', 'five_prime_utr',
                         'intron', 'gene']

        if df.shape[0] == 1:
            return df

        for feature in feature_order:
            _df = df[df['Feature'] == feature]

            if _df.shape[0] == 0:
                continue

            # feature type is not 3'UTR just get first overlap
            if feature != 'three_prime_utr':
                return _df.drop_duplicates(subset='gene_id')

            # if feature type is 3'UTR get feature with closest end
            _df['dist'] = (_df['polyA_site'] - _df['canonical_site']).abs()

            _df = _df.groupby('gene_id').apply(
                lambda xdf: xdf.loc[xdf['dist'].idxmin()])

            del _df['dist']

            return _df

    # def overlap_gene(self, gr):
    #     df_gene = self.gr[self.gr.Feature == 'gene'] \
    #                   .df.drop_duplicates(subset=['gene_id'])

    #     df_gene = df_gene[df_gene['gene_type'].isin(
    #         ['protein_coding', 'lncRNA'])]

    #     df_gene = df_gene[[*_core_columns, 'gene_id', 'gene_name']]
    #     df_gene = df_gene[df_gene['Strand'].isin(gr.Strand.cat.categories)]
    #     gr_gene = pr.PyRanges(df_gene, int64=True)

    #     gr_gene = pr.PyRanges(df_gene, int64=True)

    #     gr = gr.join(gr_gene, how='left')
    #     return gr.drop(['Start_b', 'End_b', 'Strand_b'])

    # def overlap_tes(self, gr):
    #     gr_tes = self.gr.features.tes()

    #     df_tes = gr_tes.df[_core_columns].drop_duplicates()
    #     df_tes['canonical'] = True

    #     df_tes = df_tes[df_tes['Strand'].isin(gr.Strand.cat.categories)]
    #     gr_tes = pr.PyRanges(df_tes, int64=True)

    #     return gr.join(gr_tes, how='left').drop(['End_b', 'Strand_b'])
