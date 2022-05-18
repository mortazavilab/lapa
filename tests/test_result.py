import pandas as pd
from lapa import lapa
from lapa.utils.io import read_talon_read_annot
from lapa.result import LapaResult
from conftest import fasta, gtf, chrom_sizes, read_annot


def test_LapaResult_fisher_exact_test(tmp_path):
    output_dir = tmp_path / 'lapa'

    df_read_annot = pd.read_csv(read_annot, sep='\t')
    df_read_annot2 = df_read_annot.copy()
    df_read_annot2['dataset'] += '2'

    sample_read_annot = str(tmp_path / 'sample_read_annot.tsv')

    pd.concat([df_read_annot, df_read_annot2]).to_csv(
        sample_read_annot, sep='\t', index=False)

    lapa(sample_read_annot, fasta, gtf, chrom_sizes, output_dir)

    result = LapaResult(output_dir)

    groups = {
        'hegp': ['hepg2', 'hepg22'],
        'gm12': ['gm12878', 'gm128782']
    }

    df = result.fisher_exact_test(groups)

    df.columns == [
        'odds_ratio', 'pval', 'delta_usage', 'gene_id', 'pval_adj'
    ]

    assert df[(df['pval_adj'] < 0.05) & (
        df['delta_usage'].abs() > 0.3)].shape[0]


# def test_LapaResult__count_per_group():
#     df = LapaResult._counts_per_groups(
#         pd.DataFrame({
#             'a': [1, 2, 3],
#             'b': [3, 4, 5],
#             'x': [0, 0, 1],
#             'y': [1, 1, 2]

#         }),
#         groups={
#             'c': ['a', 'b'],
#             'z': ['x', 'y']
#         }
#     )

#     pd.testing.assert_frame_equal(
#         df,
#         pd.DataFrame({
#             'c': [4, 6, 8],
#             'z': [1, 1, 3]
#         })
#     )
