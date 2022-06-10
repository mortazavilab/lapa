import logging
import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm





def _tqdm_pandas_gr():
    '''
    Adaptor for tqdm to integrate to logging
    '''
    logger = logging.getLogger('progress')
    file = logger.handlers[0].stream if logger.handlers else None

    tqdm.pandas(mininterval=5, file=file,
                bar_format='- {n_fmt} cluster is annotated...\n')


class _GenomicRegions:

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

        _core_columns = ['Chromosome', 'Start', 'End', 'Strand']
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

        df['Feature'] = df['Feature'].replace('-1', 'intergenic')
        df['gene_id'] = df['gene_id'].replace('-1', '')
        df['gene_name'] = df['gene_name'].replace('-1', '')

        _tqdm_pandas_gr()

        # if gene_id defined overlap
        _core_columns = ['Chromosome', 'Start', 'End', 'Strand']        
        df = df.groupby(_core_columns, observed=True) \
               .progress_apply(self._agg_annotation_gene) \
               .reset_index(drop=True)

        df = self.intergenic_genes(df)
        return df

    def intergenic_genes(self, df):
        '''
        Assign `gene_id` and `gene_name` to intergenic clusters.
        '''
        genes = pd.Series(range(sum(df['Feature'] == 'intergenic'))) \
                  .astype('str').add_prefix('intergenic_').index

        df.loc[df['Feature'] == 'intergenic', 'gene_id'] = genes.tolist()
        df.loc[df['Feature'] == 'intergenic', 'gene_name'] = genes.tolist()

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

            # feature type is not annotated_feature_site (UTR)
            # just get first overlap
            if feature not in {self.annotated_feature_site, 'exon'}:
                return _df.iloc[[0]]

            return _df.drop_duplicates(subset='gene_id')


class PolyAGenomicRegions(_GenomicRegions):
    '''
    Annotate polyA sites based on the genomics features.

    Args:
      gtf_file: Annotation file overlap against.

    Examples:
      Annotation of poly(A) in pyranges format:
      
      >>> regions = PolyAGenomicRegions('hg38.gtf')
      >>> gr
      +--------------+-----------+-----------+--------------+--------------+
      | Chromosome   |     Start |       End | Strand       |   polyA_site |
      | (category)   |   (int32) |   (int32) | (category)   |      (int64) |
      |--------------+-----------+-----------+--------------+--------------|
      | chr17        |   4303000 |   4303100 | -            |      4303050 |
      | chr17        |  43044826 |  43045289 | -            |     43045057 |
      | chr17        |  43046541 |  43046997 | -            |     43046769 |
      | chr17        |  43115728 |  43115767 | -            |     43057094 |
      | chr17        |  43093458 |  43093573 | -            |     43093515 |
      +--------------+-----------+-----------+--------------+--------------+
      >>> regions.annotate(gr)
      +--------------+-----------+-----------+--------------+--------------+-----------------+--------------------+--------------+------------------+
      | Chromosome   |     Start |       End | Strand       |   polyA_site | Feature         | gene_id            | gene_name    |   annotated_site |
      | (category)   |   (int32) |   (int32) | (category)   |      (int64) | (object)        | (object)           | (object)     |          (int64) |
      |--------------+-----------+-----------+--------------+--------------+-----------------+--------------------+--------------+------------------|
      | chr17        |   4303000 |   4303100 | -            |      4303050 | intergenic      | intergenic_0       | intergenic_0 |               -1 |
      | chr17        |  43044826 |  43045289 | -            |     43045057 | three_prime_utr | ENSG00000012048.23 | BRCA1        |         43044294 |
      | chr17        |  43046541 |  43046997 | -            |     43046769 | intron          | ENSG00000012048.23 | BRCA1        |               -1 |
      | chr17        |  43093458 |  43093573 | -            |     43093515 | three_prime_utr | ENSG00000012048.23 | BRCA1        |         43091434 |
      | chr17        |  43115728 |  43115767 | -            |     43057094 | exon            | ENSG00000012048.23 | BRCA1        |               -1 |
      +--------------+-----------+-----------+--------------+--------------+-----------------+--------------------+--------------+------------------+
    '''

    def __init__(self, gtf_file):
        feature_order = ['three_prime_utr', 'exon',
                         'intron', 'five_prime_utr', 'gene']
        super().__init__(gtf_file, feature_order, 'three_prime_utr')

    def annotated_site(self, df):
        return np.where(
            df['Feature'] == self.annotated_feature_site,
            np.where(df['Strand'] == '-', df['Start'], df['End']), -1)


class TssGenomicRegions(_GenomicRegions):
 
    '''
    Annotate tss sites based on the genomics features.

    Args:
      gtf_file: Annotation file overlap against.

    Examples:
      Annotation of tss in pyranges format:
      
      >>> regions = TssGenomicRegions('hg38.gtf')
      >>> gr
      +--------------+-----------+-----------+--------------+--------------+
      | Chromosome   |     Start |       End | Strand       |   polyA_site |
      | (category)   |   (int32) |   (int32) | (category)   |      (int64) |
      |--------------+-----------+-----------+--------------+--------------|
      | chr17        |   4303000 |   4303100 | -            |      4303050 |
      | chr17        |  43044826 |  43045289 | -            |     43045057 |
      | chr17        |  43046541 |  43046997 | -            |     43046769 |
      | chr17        |  43115728 |  43115767 | -            |     43057094 |
      | chr17        |  43093458 |  43093573 | -            |     43093515 |
      +--------------+-----------+-----------+--------------+--------------+
      >>> regions.annotate(gr)
      +--------------+-----------+-----------+--------------+--------------+-----------------+--------------------+--------------+------------------+
      | Chromosome   |     Start |       End | Strand       |   polyA_site | Feature         | gene_id            | gene_name    |   annotated_site |
      | (category)   |   (int32) |   (int32) | (category)   |      (int64) | (object)        | (object)           | (object)     |          (int64) |
      |--------------+-----------+-----------+--------------+--------------+-----------------+--------------------+--------------+------------------|
      | chr17        |   4303000 |   4303100 | -            |      4303050 | intergenic      | intergenic_0       | intergenic_0 |               -1 |
      | chr17        |  43044826 |  43045289 | -            |     43045057 | five_prime_utr  | ENSG00000012048.23 | BRCA1        |         43044294 |
      | chr17        |  43046541 |  43046997 | -            |     43046769 | intron          | ENSG00000012048.23 | BRCA1        |               -1 |
      | chr17        |  43093458 |  43093573 | -            |     43093515 | five_prime_utr  | ENSG00000012048.23 | BRCA1        |         43091434 |
      | chr17        |  43115728 |  43115767 | -            |     43057094 | exon            | ENSG00000012048.23 | BRCA1        |               -1 |
      +--------------+-----------+-----------+--------------+--------------+-----------------+--------------------+--------------+------------------+
    '''

    def __init__(self, gtf_file):
        feature_order = ['five_prime_utr', 'exon',
                         'intron', 'three_prime_utr', 'gene']
        super().__init__(gtf_file, feature_order, 'five_prime_utr')

    def annotated_site(self, df):
        return np.where(
            df['Feature'] == self.annotated_feature_site,
            np.where(df['Strand'] == '-', df['End'], df['Start']), -1)
