import pandas as pd


polyA_signal_seqs = [
    'AATAAA', 'ATTAAA', 'TATAAA', 'AGTAAA',
    'AATACA', 'CATAAA', 'AATATA', 'GATAAA',
    'AATGAA', 'AATAAT', 'AAGAAA', 'ACTAAA',
    'AATAGA', 'ATTACA', 'AACAAA', 'ATTATA',
    'AACAAG', 'AATAAG'
]
'''
List[str]: the list of poly(A) signals
that are annotated for single 3' end processing sites
in the region -60 to +10 nt around them
According to Gruber et al., 2016, Genome Research.
'''


def pad_series(series, pad_size=2):
    '''
    Pad pandas series with zeros at the start ande end.

    Args:
      series: pd.Series
      pad_size: 0
    '''
    return pd.concat([
        pd.Series([0] * pad_size),
        pd.Series(series),
        pd.Series([0] * pad_size)
    ])
