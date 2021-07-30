import pandas as pd
import pyranges as pr


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
        df = df[['Chromosome', 'Start', 'End', 'Strand', 'Feature']]
        return pr.PyRanges(df, int64=True)

    def annotate(self, gr, features=None, single=False):
        gr_ann = pr.PyRanges(gr.df, int64=True).join(
            self.features(features),
            strandedness='same', how='left')
        df = gr_ann.df
        # need to avoid duplication of groupby for category
        df['Chromosome'] = df['Chromosome'].astype(str)
        df['Strand'] = df['Strand'].astype(str)
        df['Feature'] = df['Feature'].replace('-1', 'intergenic')

        df = df[['Chromosome', 'Start', 'End', 'Strand', 'Feature']] \
            .groupby(['Chromosome', 'Start', 'End', 'Strand']) \

        df = df.agg(set)
        if single:
            order = ['three_prime_utr', 'five_prime_utr', 'UTR',
                     'exon', 'intron', 'gene', 'intergenic']

            def _get_feature(x):
                for i in order:
                    if i in x:
                        return i
            df['Feature'] = df['Feature'].map(_get_feature)
        else:
            df['Feature'] = df['Feature'].map(lambda x: ','.join(sorted(x)))

        return pr.PyRanges(gr.df.set_index(
            ['Chromosome', 'Start', 'End', 'Strand']).join(df).reset_index(), int64=True)
