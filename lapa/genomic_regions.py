import logging
import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm


_core_columns = ['Chromosome', 'Start', 'End', 'Strand']


def _tqdm_pandas_gr():
    '''
    Adaptor for tqdm to integrate to logging
    '''
    logger = logging.getLogger('progress')
    file = logger.handlers[0].stream if logger.handlers else None

    tqdm.pandas(mininterval=5, file=file,
                bar_format='- {n_fmt} cluster is annotated...\n')


class GenomicRegions:

    def __init__(self, gtf_file, feature_order, annotated_feature_site):
        if isinstance(gtf_file, str):
            self.gr = pr.read_gtf(gtf_file)
        elif isinstance(gtf_file, pr.PyRanges):
            self.gr = gtf_file
        else:
            raise ValueError('invalid gtf_file argument')

        assert annotated_feature_site in self.gr.Feature.unique(), \
            f'{annotated_feature_site} not in the gtf file. ' \
            'If you gtf file dont have five_prime_utr, three_prime_utr' \
            ' but just UTR run https://github.com/MuhammedHasan/gencode_utr_fix'

        self.feature_order = feature_order
        self.annotated_feature_site = annotated_feature_site

    def annotated_site(self, df):
        raise NotImplementedError()

    def features(self, features: set = None) -> pr.PyRanges:
        features = features or set(self.feature_order)

        gr = self.gr[self.gr.Feature.isin(features)]
        df = gr.df

        if 'intron' in features:
            gr_introns = self.gr.features.introns()
            df = pd.concat([df, gr_introns.df])

        df = df[[*_core_columns, 'Feature', 'gene_id', 'gene_name']]

        df['annotated_site'] = self.annotated_site(df)

        return df

    def annotate(self, gr, features=None):

        df_gtf = self.features(features)
        df_gtf = df_gtf[df_gtf['Strand'].isin(gr.Strand.cat.categories)]
        gr_gtf = pr.PyRanges(df_gtf, int64=True)

        gr_ann = pr.PyRanges(gr.df, int64=True).join(
            gr_gtf, strandedness='same', how='left') \
            .drop(['Start_b', 'End_b', 'Strand_b'])
        df = gr_ann.df.drop_duplicates()

        # need to avoid duplication of groupby for category
        # df['Chromosome'] = df['Chromosome'].astype(str)
        # df['Strand'] = df['Strand'].astype(str)

        df['Feature'] = df['Feature'].replace('-1', 'intergenic')
        df['gene_id'] = df['gene_id'].replace('-1', '')
        df['gene_name'] = df['gene_name'].replace('-1', '')

        _tqdm_pandas_gr()

        # if gene_id defined overlap
        df = df.groupby(_core_columns, observed=True) \
               .progress_apply(self._agg_annotation_gene) \
               .reset_index(drop=True)

        return df

    def _agg_annotation_gene(self, df):
        feature_order = self.feature_order

        if df.shape[0] == 1:
            return df

        for feature in feature_order:
            _df = df[df['Feature'] == feature]

            if _df.shape[0] == 0:
                continue

            if _df.shape[0] == 1:
                return _df

            # feature type is not 3'UTR just get first overlap
            if feature not in {'three_prime_utr', 'exon'}:
                return _df.iloc[[0]]

            return _df.drop_duplicates(subset='gene_id')


class PolyAGenomicRegions(GenomicRegions):

    def __init__(self, gtf_file):
        feature_order = ['three_prime_utr', 'exon',
                         'intron', 'five_prime_utr', 'gene']
        super().__init__(gtf_file, feature_order, 'three_prime_utr')

    def annotated_site(self, df):
        return np.where(
            df['Feature'] == self.annotated_feature_site,
            np.where(df['Strand'] == '-', df['Start'], df['End']), -1)


class TssGenomicRegions(GenomicRegions):

    def __init__(self, gtf_file):
        feature_order = ['five_prime_utr', 'exon',
                         'intron', 'three_prime_utr', 'gene']
        super().__init__(gtf_file, feature_order, 'five_prime_utr')

    def annotated_site(self, df):
        return np.where(
            df['Feature'] == self.annotated_feature_site,
            np.where(df['Strand'] == '-', df['End'], df['Start']), -1)
