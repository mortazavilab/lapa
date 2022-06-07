from pathlib import Path
import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from lapa.utils.io import read_polyA_cluster, read_tss_cluster
from lapa.replication import replication_rate

_core_cols = ['Chromosome', 'Start', 'End', 'Strand']


class _LapaResult:

    def __init__(self, path, replicated=True, prefix=''):
        self.lapa_dir = Path(path)
        self.replicated = replicated
        self.prefix = prefix

        self._samples = None
        self._datasets = None

    def _read_cluster(self, path):
        raise NotImplementedError()

    @property
    def samples(self):
        if self._samples is None:
            self._samples = [
                i.stem.replace('.bed', '')
                for i in self.sample_dir.iterdir()
            ]
        return self._samples

    @property
    def datasets(self):
        if self._datasets is None:
            self._datasets = [
                i.stem.replace('.bed', '')
                for i in self.dataset_dir.iterdir()
            ]
        return self._datasets

    @property
    def count_dir(self):
        return self.lapa_dir / 'counts'

    @property
    def sample_dir(self):
        if not self.replicated:
            return self.lapa_dir / 'raw_sample'
        return self.lapa_dir / 'sample'

    @property
    def dataset_dir(self):
        return self.lapa_dir / 'dataset'

    @property
    def cluster_path(self):
        if not self.replicated:
            return self.lapa_dir / (f'raw_{self.prefix}_clusters.bed')
        return self.lapa_dir / (f'{self.prefix}_clusters.bed')

    @property
    def dataset_path(self):
        if not self.replicated:
            return self.lapa_dir / (f'raw_{self.prefix}_clusters.bed')
        return self.lapa_dir / (f'{self.prefix}_clusters.bed')

    def read_clusters(self, filter_intergenic=True):
        df = self._read_cluster(self.cluster_path)

        if filter_intergenic:
            df = df[df['Feature'] != 'intergenic']

        return self._set_index(df)

    def read_sample(self, sample, filter_intergenic=True):
        if sample not in self.samples:
            raise ValueError(
                'sample `%s` does not exist in directory' % sample)

        df = self._read_cluster(self.sample_dir / ('%s.bed' % sample))

        if filter_intergenic:
            df = df[df['Feature'] != 'intergenic']

        return self._set_index(df)

    def read_dataset(self, dataset):
        df = self._read_cluster(self.dataset_dir / ('%s.bed' % dataset))
        return self._set_index(df)

    def _set_index(self, df):
        df['name'] = df['Chromosome'] + ':' \
            + df[f'{self.prefix}_site'].astype('str') + ':' \
            + df['Strand'].astype('str')
        return df.set_index('name')

    def read_counts(self, sample=None, strand=None):
        sample = sample or 'all'

        if strand == '+':
            df = pr.read_bigwig(
                str(self.count_dir / (f'{sample}_{self.prefix}_counts_pos.bw'))
            ).df
            df['Strand'] = '+'
            return df.rename(columns={'Value': 'count'})
        elif strand == '-':
            df = pr.read_bigwig(
                str(self.count_dir / (f'{sample}_{self.prefix}_counts_neg.bw'))
            ).df
            df['Strand'] = '-'
            return df.rename(columns={'Value': 'count'})
        else:
            return pd.concat([
                self.read_counts(sample, strand='+'),
                self.read_counts(sample, strand='-')
            ])

    def attribute(self, field, samples=None):
        '''
        Read attribute of samples as dataframe
        '''
        samples = samples or self.samples

        df = pd.concat([
            self.read_sample(sample)
            .drop_duplicates(_core_cols)
            .rename(columns={field: sample})[sample]
            for sample in self.samples
        ], axis=1).sort_index()
        df.index = df.index.rename('site')
        return df

    def counts(self, samples=None):
        return self.attribute('count', samples=samples)

    def total_counts(self, samples=None):
        return self.counts(samples=samples).sum(axis=1)

    def gene_id(self):
        return self.attribute('gene_id').apply(
            lambda row: row[~row.isna()][0], axis=1)

    @staticmethod
    def _agg_per_groups(df, groups, agg_func):
        return pd.DataFrame({
            k: df[v].agg(agg_func, axis=1)
            for k, v in groups.items()
        })

    def _k_n(self, groups, min_gene_count):
        counts = self.counts()

        k = self._agg_per_groups(counts, groups, 'sum')
        n = self.attribute('gene_count')
        n = self._agg_per_groups(n, groups, 'sum')

        filter_rows = ((n > 0).sum(axis=1) > 1) \
            & (k != n).any(axis=1) \
            & (n > min_gene_count).any(axis=1)

        return k[filter_rows], n[filter_rows]

    def replication_rate(self, samples=None, num_samples=2, min_score=1):
        '''
        Calculate replication rate of samples
        '''
        samples = samples or self.samples

        return replication_rate({
            i: self.read_sample(i, filter_intergenic=False)
            for i in samples
        }, score_column='count', num_samples=num_samples,
            min_score=min_score)

    def plot_replication_rate(self, samples=None, num_samples=2,
                              min_score=1, line_kws=None):
        import seaborn as sns

        samples = samples or self.samples
        line_kws = line_kws or dict()

        df = self.replication_rate(samples, num_samples=num_samples,
                                   min_score=min_score)
        df['rank'] = df['score'].rank(ascending=False)
        return sns.lineplot(data=df, x='rank', y='replication', **line_kws)

    def fisher_exact_test(self, groups, min_gene_count=10,
                          correction_method='fdr_bh'):
        '''
        Fisher-exact test for sites.

        Args:
            groups (Dict[str, List[str]]): dict of two elements
                as assinging groups. Two keys are group names
                and values are list of keys annotating samples
                belong to each group.
            min_gene_count (int): Number of reads in the gene to
                be consider in analysis.
            correction_method (str): multiple testing correction method.
                methods in `statsmodels.stats.multitest.multipletests`
                are valid.
        '''
        if len(groups) != 2:
            raise 'Two groups are need for fisher_exact test'

        k, n = self._k_n(groups, min_gene_count)

        odds, pvals = zip(*[
            fisher_exact([
                [_k[0],         _k[1]],
                [_n[0] - _k[0], _n[1] - _k[1]]
            ])
            for _k, _n in tqdm(zip(k.values, n.values), total=k.shape[0])
        ])

        sites = k.index.tolist()

        usage_dif = self._agg_per_groups(
            self.attribute('usage').loc[sites], groups, 'mean')
        groups = usage_dif.columns
        usage_dif = usage_dif[groups[0]] - usage_dif[groups[1]]

        df = pd.DataFrame({
            'odds_ratio': odds,
            'pval': pvals,
            'delta_usage': usage_dif,
            'gene_id': self.gene_id().loc[sites]
        }, index=k.index.tolist())

        df['pval_adj'] = multipletests(df['pval'], method=correction_method)[1]
        return df

    # beta-binomial test

    def beta_binomial_test(self, min_gene_count=10, theta=0.001, max_iter=1000):
        '''
        P-values based on betabinomial test.
        '''
        from betabinomial import BetaBinomial, pval_adj

        # merge replicates into one count
        counts = self.counts()
        k = counts.values
        n = self.attribute('gene_count').values

        filter_rows = ((n > 0).sum(axis=1) > 1) \
            & (k != n).any(axis=1) \
            & (n > min_gene_count).any(axis=1)

        n = n[filter_rows]
        k = k[filter_rows]
        bb = BetaBinomial().infer(k, n, theta=theta, max_iter=max_iter)

        # recalculate usage k/n
        usage = self.attribute('usage')
        usage = usage[filter_rows]

        gene_id = self.gene_id()
        gene_id = gene_id[filter_rows]
        gene_id = np.repeat(gene_id.values.reshape((-1, 1)), 4, axis=1)

        sites = counts.index
        sites = sites[filter_rows]

        cols = {
            'usage': usage,
            'delta_usage': usage - bb.beta_mean(),
            'count': k,
            'expected_count': bb.mean(n),
            'gene_count': n,
            'gene_id': gene_id,
            'pval': bb.pval(k, n),
            'z_score': bb.z_score(k, n),
            'logfc': bb.log_fc(k, n)
        }
        cols['padj'] = pval_adj(np.nan_to_num(cols['pval'], nan=1))

        df = pd.concat([
            pd.DataFrame(v, index=sites, columns=self.samples)
            .reset_index()
            .melt(id_vars='polya_site', var_name='sample', value_name=col)
            .set_index(['polya_site', 'sample'])
            for col, v in cols.items()
        ], axis=1)

        return df


class LapaResult(_LapaResult):

    def __init__(self, path, replicated=True):
        super().__init__(path, replicated, 'polyA')

    def _read_cluster(self, path):
        return read_polyA_cluster(path)

    def read_clusters(self, filter_intergenic=True,
                      filter_internal_priming=True):
        df = super().read_clusters(self.cluster_path)

        if filter_internal_priming:
            df = df[~(
                (df['fracA'] > 7) &
                (df['signal'] == 'None@None')
            )]

        return df


class LapaTssResult(_LapaResult):

    def __init__(self, path, replicated=True):
        super().__init__(path, replicated, 'tss')

    def _read_cluster(self, path):
        return read_tss_cluster(path)
