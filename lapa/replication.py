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


def replication_rate(samples, score_column='Score', rolling_size=1000):
    '''
    Returns replication rate for samples for score_column.

    Args:
      sample: Dictionary of sample name as key and dataframe with
        columns of 'Chromosome', 'Start', 'End', 'Strand' as value
      score_column: Column to calculate rank samples
      rolling_size: Rolling size in replication rate calculation
    '''
    if len(samples) < 2:
        raise ValueError('Cannot calcultate replication rate '
                         'because number of samples are smaller than 2')

    df = agg_sample_cols(samples, score_column)

    df = pd.DataFrame({
        'replication': (df > 0).sum(axis=1) / df.shape[1],
        'score': df.sum(axis=1)
    }).sort_values('score', ascending=False).rolling(rolling_size).mean()

    return df


def replication_dataset(samples, score_column='Score', rolling_size=1000,
                        min_replication_rate=0.75):
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

    df = replication_rate(samples, score_column, rolling_size)
    indexes = df[df['replication'].fillna(1) >= min_replication_rate].index

    return {
        sample: _df[_df.set_index(_core_cols).index.isin(indexes)]
        for sample, _df in samples.items()
    }
