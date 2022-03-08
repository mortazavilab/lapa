from pathlib import Path
import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
# from betabinomial import BetaBinomial, pval_adj
from lapa.utils.io import read_apa_sample, read_polyA_cluster


_core_cols = ['Chromosome', 'Start', 'End', 'Strand']


class LapaResult:

    def __init__(self, path, tpm_cutoff=1):
        self.lapa_dir = Path(path)
        self.tpm_cutoff = tpm_cutoff
        self.samples = [
            i.stem.replace('_apa', '')
            for i in self.lapa_dir.iterdir()
            if i.stem.endswith('_apa')
        ]

    def read_apa(self, sample):
        if sample not in self.samples:
            raise ValueError(
                'sample `%s` does not exist in directory' % sample)
        df = read_apa_sample(self.lapa_dir / ('%s_apa.bed' % sample))
        return self._filter_tpm(self._set_index(df))

    def read_cluster(self, filter_internal_priming=True):
        df = read_polyA_cluster(self.lapa_dir / 'polyA_clusters.bed')
        if filter_internal_priming:
            df = df[(
                ~(
                    (df['fracA'] > 7) &
                    (df['signal'] == 'None@None')
                )) | (df['canonical_site'] != -1)
            ]
        return self._filter_tpm(self._set_index(df))

    def read_counts(self, sample=None, strand=None):
        sample = sample or 'all'

        if strand == '+':
            df = pr.read_bigwig(
                str(self.lapa_dir / ('%s_tes_counts_pos.bw' % sample))).df
            df['Strand'] = '+'
            return df.rename(columns={'Value': 'count'})
        elif strand == '-':
            df = pr.read_bigwig(
                str(self.lapa_dir / ('%s_tes_counts_neg.bw' % sample))).df
            df['Strand'] = '-'
            return df.rename(columns={'Value': 'count'})
        else:
            return pd.concat([
                self.read_counts(sample, strand='+'),
                self.read_counts(sample, strand='-')
            ])

    def attribute(self, field):
        df = pd.concat([
            self.read_apa(sample)
            .drop_duplicates(_core_cols)
            .rename(columns={field: sample})[sample]
            for sample in self.samples
        ], axis=1).sort_index()
        df.index = df.index.rename('polya_site')
        return df

    def counts(self):
        return self.attribute('count')

    def total_counts(self):
        return self.counts().sum(axis=1)

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

    def fisher_exact_test(self, groups, min_gene_count=10,
                          correction_method='fdr_bh'):
        '''
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

        polya_sites = k.index.tolist()

        usage_dif = self._agg_per_groups(
            self.attribute('usage').loc[polya_sites], groups, 'mean')
        groups = usage_dif.columns
        usage_dif = usage_dif[groups[0]] - usage_dif[groups[1]]

        df = pd.DataFrame({
            'odds_ratio': odds,
            'pval': pvals,
            'delta_usage': usage_dif,
            'gene_id': self.gene_id().loc[polya_sites]
        }, index=k.index.tolist())

        df['pval_adj'] = multipletests(df['pval'], method=correction_method)[1]
        return df

    # beta-binomial test
    # def stats_testing(self, groups=None, min_gene_count=10,
    #                   theta=0.001, max_iter=1000):
    #     '''
    #     P-values based on betabinomial test.
    #     '''
    #     # merge replicates into one count
    #     counts = self.counts()
    #     k = counts.values
    #     n = self.attribute('gene_count').values

    #     filter_rows = ((n > 0).sum(axis=1) > 1) \
    #         & (k != n).any(axis=1) \
    #         & (n > min_gene_count).any(axis=1)

    #     n = n[filter_rows]
    #     k = k[filter_rows]
    #     bb = BetaBinomial().infer(k, n, theta=theta, max_iter=max_iter)

    #     # recalculate usage k/n
    #     usage = self.attribute('usage')
    #     usage = usage[filter_rows]

    #     gene_id = self.gene_id()
    #     gene_id = gene_id[filter_rows]
    #     gene_id = np.repeat(gene_id.values.reshape((-1, 1)), 4, axis=1)

    #     sites = counts.index
    #     sites = sites[filter_rows]ll

    #     cols = {
    #         'usage': usage,
    #         'delta_usage': usage - bb.beta_mean(),
    #         'count': k,
    #         'expected_count': bb.mean(n),
    #         'gene_count': n,
    #         'gene_id': gene_id,
    #         'pval': bb.pval(k, n),
    #         'z_score': bb.z_score(k, n),
    #         'logfc': bb.log_fc(k, n)
    #     }
    #     cols['padj'] = pval_adj(np.nan_to_num(cols['pval'], nan=1))

    #     df = pd.concat([
    #         pd.DataFrame(v, index=sites, columns=self.samples)
    #         .reset_index()
    #         .melt(id_vars='polya_site', var_name='sample', value_name=col)
    #         .set_index(['polya_site', 'sample'])
    #         for col, v in cols.items()
    #     ], axis=1)

    #     return df

    @staticmethod
    def _set_index(df):
        df['name'] = df['Chromosome'] + ':' \
            + df['polyA_site'].astype('str') + ':' \
            + df['Strand'].astype('str')
        return df.set_index('name')

    def _filter_tpm(self, df):
        return df[df['tpm'] >= self.tpm_cutoff]
