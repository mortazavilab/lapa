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
        df = df[[*_core_columns, 'Feature']]
        return df

    def annotate(self, gr, features=None, single=False):
        df_gtf = self.features(features)
        df_gtf = df_gtf[df_gtf['Strand'].isin(gr.Strand.cat.categories)]
        gr_gtf = pr.PyRanges(df_gtf, int64=True)

        gr_ann = pr.PyRanges(gr.df, int64=True).join(
            gr_gtf, strandedness='same', how='left')
        df = gr_ann.df
        # need to avoid duplication of groupby for category
        df['Chromosome'] = df['Chromosome'].astype(str)
        df['Strand'] = df['Strand'].astype(str)

        # df['Feature'] = df['Feature'].replace('-1', 'intergenic')

        df = df[[*_core_columns, 'Feature']] \
            .groupby(_core_columns).agg(set)

        if single:
            order = ['three_prime_utr', 'five_prime_utr', 'UTR',
                     'exon', 'intron', 'gene']

            def _get_feature(x):
                for i in order:
                    if i in x:
                        return i
            df['Feature'] = df['Feature'].map(_get_feature)
        else:
            df['Feature'] = df['Feature'].map(lambda x: ','.join(sorted(x)))

        return pr.PyRanges(gr.df.set_index(
            _core_columns).join(df).reset_index(), int64=True)

    def overlap_gene(self, gr):
        df_gene = self.gr[self.gr.Feature == 'gene'].df

        df_gene = df_gene[df_gene['gene_type'].isin(
            ['protein_coding', 'lncRNA'])]

        df_gene = df_gene[[*_core_columns, 'gene_id', 'gene_name']]
        df_gene = df_gene[df_gene['Strand'].isin(gr.Strand.cat.categories)]
        gr_gene = pr.PyRanges(df_gene, int64=True)

        gr_gene = pr.PyRanges(df_gene, int64=True)
        gr = gr.join(gr_gene, how='left')

        return gr.drop(['Start_b', 'End_b', 'Strand_b'])

    def overlap_tes(self, gr):
        gr_tes = self.gr.features.tes()

        df_tes = gr_tes.df[_core_columns].drop_duplicates()
        df_tes['canonical'] = True

        df_tes = df_tes[df_tes['Strand'].isin(gr.Strand.cat.categories)]
        gr_tes = pr.PyRanges(df_tes, int64=True)

        return gr.join(gr_tes, how='left').drop(['End_b', 'Strand_b'])
