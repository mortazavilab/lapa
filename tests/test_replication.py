import pytest
import numpy as np
import pandas as pd
from lapa.replication import agg_sample_cols, replication_rate, \
    replication_dataset


@pytest.fixture
def samples_replicates():
    return {
        'a': pd.DataFrame({
            'Chromosome': ['chr1', 'chr2', 'chr1'],
            'Start': [10, 1000, 20],
            'End': [100, 1010, 120],
            'Strand': ['+', '-', '+'],
            'Score': [10, 2, 1]
        }),
        'b': pd.DataFrame({
            'Chromosome': ['chr1', 'chr2', 'chr1'],
            'Start': [10, 1000, 100],
            'End': [100, 1010, 150],
            'Strand': ['+', '-', '-'],
            'Score': [20, 2, 1]
        })
    }


def test_agg_sample_cols(samples_replicates):

    df = agg_sample_cols(samples_replicates)

    pd.testing.assert_frame_equal(
        df.reset_index(),
        pd.DataFrame({
            'Chromosome': ['chr1', 'chr1', 'chr1', 'chr2'],
            'Start': [10, 20, 100, 1000],
            'End': [100, 120, 150, 1010],
            'Strand': ['+', '+', '-', '-'],
            'a': [10.0, 1.0, 0.0, 2.0],
            'b': [20.0, 0.0, 1.0, 2.0]
        })
    )


def test_replication_rate(samples_replicates):
    df = replication_rate(samples_replicates, rolling_size=2)

    pd.testing.assert_frame_equal(
        df.reset_index(),
        pd.DataFrame({
            'Chromosome': ['chr1', 'chr2', 'chr1', 'chr1'],
            'Start': [10, 1000, 20, 100],
            'End': [100, 1010, 120, 150],
            'Strand': ['+', '-', '+', '-'],
            'replication': [np.nan, 1.0, 0.5, 0],
            'score': [np.nan, 17.0, 2.5, 1.0]
        })
    )


def test_replication_dataset(samples_replicates):

    samples = replication_dataset(samples_replicates, rolling_size=2)

    pd.testing.assert_frame_equal(
        samples['a'],
        pd.DataFrame({
            'Chromosome': ['chr1', 'chr2'],
            'Start': [10, 1000],
            'End': [100, 1010],
            'Strand': ['+', '-'],
            'Score': [10, 2]
        })
    )

    pd.testing.assert_frame_equal(
        samples['b'],
        pd.DataFrame({
            'Chromosome': ['chr1', 'chr2'],
            'Start': [10, 1000],
            'End': [100, 1010],
            'Strand': ['+', '-'],
            'Score': [20, 2]
        })
    )
