import pandas as pd


chroms_chr = {
    'chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
    'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr3', 'chr4',
    'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX', 'chrY'
}

chroms = {
    '1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
    '2', '3', '4', '5', '6', '7', '8', '9', 'X', 'Y'
}


polyA_signal_seqs = [
    'AATAAA', 'ATTAAA', 'TATAAA', 'AGTAAA',
    'AATACA', 'CATAAA', 'AATATA', 'GATAAA',
    'AATGAA', 'AATAAT', 'AAGAAA', 'ACTAAA',
    'AATAGA', 'ATTACA', 'AACAAA', 'ATTATA',
    'AACAAG', 'AATAAG'
]


def pad_series(series, pad_size=2):
    return pd.concat([
        pd.Series([0] * pad_size),
        pd.Series(series),
        pd.Series([0] * pad_size)
    ])


def filter_main_chroms(df):
    return df[df['Chromosome'].isin(chroms)]
