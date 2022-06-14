import pandas as pd
import pyranges as pr
from lapa.utils.io import bw_from_pyranges


def test_bw_from_pyranges(tmp_path):
    
    gr = pr.PyRanges(pd.DataFrame({
        'Chromosome': ['chrM', 'chrM', 'chr1', 'chr1'],
        'Start': [-1, 0, 138, 20],
        'End': [0, 1, 149, 21],
        'Strand': ["-", '-', "+", '-'],
        'count': [1, 2, 3, 10]
    }))

    bw_from_pyranges(
        gr, 'count',
        chrom_sizes={'chr1': 1000, 'chr2': 10, 'chrM': 10},
        bw_pos_file = str(tmp_path / 'pos.bw'),
        bw_neg_file = str(tmp_path / 'neg.bw')
    )
