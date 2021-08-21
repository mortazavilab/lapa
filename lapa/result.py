from pathlib import Path
import pandas as pd
import pyranges as pr
from lapa.utils.io import read_apa_sample, read_polyA_cluster


class LapaResult:

    def __init__(self, path):
        self.lapa_dir = Path(path)
        self.samples = [
            i.stem.replace('_apa', '')
            for i in self.lapa_dir.iterdir()
            if i.stem.endswith('_apa')
        ]

    def read_apa(self, sample):
        if sample not in self.samples:
            raise ValueError(
                'sample `%s` does not exist in directory' % sample)
        return read_apa_sample(self.lapa_dir / ('%s_apa.bed' % sample))

    def read_cluster(self):
        return read_polyA_cluster(self.lapa_dir / 'polyA_clusters.bed')

    def read_counts(self, sample=None, strand=None):
        sample = sample or 'all'

        if strand == '+':
            return pr.read_bigwig(
                str(self.lapa_dir / ('%s_tes_counts_pos.bw' % sample))).df
        elif strand == '-':
            return pr.read_bigwig(
                str(self.lapa_dir / ('%s_tes_counts_neg.bw' % sample))).df
        else:
            return pd.concat([
                self.read_counts(sample, strand='+'),
                self.read_counts(sample, strand='-')
            ])

        # gr_pos = pr.read_bigwig(str(lapa_dir / 'all_tes_counts_pos.bw')).df
        # gr_pos['Strand'] = '+'

        # gr = pr.PyRanges(pd.concat([gr_pos, gr_neg]))
