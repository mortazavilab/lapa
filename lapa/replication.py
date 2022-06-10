import logging
import pandas as pd


_core_cols = ['Chromosome', 'Start', 'End', 'Strand']


def agg_sample_cols(samples, column='Score'):
    '''
    Concats the column from all the dfs into one dataframe
    where column named with sample name in the dict key.

    Args:
      sample: Dictionary of sample name as key and dataframe with
        columns of 'Chromosome', 'Start', 'End', 'Strand' as value
      columns: Column to aggreate cross samples
    '''
    return pd.concat([
        df.drop_duplicates(_core_cols).set_index(_core_cols)
        .rename(columns={column: sample})[sample]
        for sample, df in samples.items()
    ], axis=1).sort_index().fillna(0)


def replication_rate(samples, score_column='Score', rolling_size=1000,
                     num_samples=2, min_score=1):
    '''
    Returns replication rate for samples for score_column.

    Args:
      sample: Dictionary of sample name as key and dataframe with
        columns of 'Chromosome', 'Start', 'End', 'Strand' as value
      score_column: Column to calculate rank samples
      rolling_size: Rolling size in replication rate calculation
      min_sample: Number of samples which region need to be
        observed to be replicated.
      min_score: Minimum score needed to recognize region as expressed.
    '''
    if len(samples) < 2:
        raise ValueError('Cannot calcultate replication rate '
                         'because number of samples are smaller than 2')

    if len(samples) < num_samples:
        logging.getLogger('warning').warning(
            f'Number of samples need for replication {num_samples}'
            ' is lower than {len(samples)}')
        num_samples = len(samples)

    df = agg_sample_cols(samples, score_column)

    df = pd.DataFrame({
        'replication': (df >= min_score).sum(axis=1) >= num_samples,
        'score': df.sum(axis=1)
    }).sort_values('score', ascending=False).rolling(rolling_size).mean()

    return df


def replication_dataset(samples, score_column='Score', rolling_size=1000,
                        min_replication_rate=0.95, num_sample=2, min_score=1):
    '''
    Calculates replication rate and filters samples based on
    given replication rate.

    Args:
      sample: Dictionary of sample name as key and dataframe with
        columns of 'Chromosome', 'Start', 'End', 'Strand' as value
      score_column: Column to calculate rank samples
      rolling_size: Rolling size in replication rate calculation
      min_replication_rate: minimum replication rate to filter samples
    '''
    if len(samples) < 2:
        raise ValueError('Cannot calcultate replication rate '
                         'because number of samples are smaller than 2')

    df = replication_rate(samples, score_column,
                          rolling_size, num_sample, min_score)
    indexes = df[df['replication'].fillna(1) >= min_replication_rate].index

    return {
        sample: _df[_df.set_index(_core_cols).index.isin(indexes)]
        for sample, _df in samples.items()
    }
