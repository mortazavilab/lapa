from pathlib import Path
import pandas as pd
import pyranges as pr
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
        return pd.concat([
            self.read_apa(sample)
            .drop_duplicates(_core_cols)
            .rename(columns={field: sample})[sample]
            for sample in self.samples
        ], axis=1).sort_index()

    def counts(self):
        return self.attribute('count')

    def total_counts(self):
        return self.counts().sum(axis=1)

    @staticmethod
    def _set_index(df):
        df['name'] = df['Chromosome'] + ':' \
            + df['polyA_site'].astype('str') + ':' \
            + df['Strand'].astype('str')
        return df.set_index('name')

    def _filter_tpm(self, df):
        return df[df['tpm'] > self.tpm_cutoff]
