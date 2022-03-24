import pytest
import pysam
import pyBigWig
import pandas as pd
import numpy as np
from lapa.count import PolyaTailCounter, FivePrimeCounter, ThreePrimeCounter, \
    TesMultiCounter, TssMultiCounter
from conftest import quantseq_gm12_bam, short_bam, chr17_chrom_sizes, read_annot


def load_reads(seq, cigar, flag, pos):
    bam = pysam.AlignmentFile(quantseq_gm12_bam)
    read = {
        'name': 'NS500169:839:HMNMLAFX2:1:11309:11029:9457',
        'flag': flag,
        'ref_name': 'chr17',
        'ref_pos': pos,
        'map_quality': '255',
        'cigar': cigar,
        'next_ref_name': '*',
        'next_ref_pos': '0',
        'length': '0',
        'seq': seq,
        'qual': '5' * len(seq),
        'tags': []
    }
    return pysam.AlignedSegment.from_dict(read, bam.header)


def test_detect_polyA():
    seq = 'CGCAATATGACGTTTCCATTTACTTTGGATTATATGTCATTATAAATATTAACAAATAAGACTTAAAAAGGACACCTTCGGGTAGGTCAGACCAAAATACAAAACTTGTCTGTGGGACTGCAGTTTGGA'
    cigar = '127M2S'
    flag = '16'
    pos = '143694'
    read = load_reads(seq, cigar, flag, pos)

    polyA_site, tail_len, percent_A = PolyaTailCounter.detect_polyA_tail(read)

    assert polyA_site is None
    assert tail_len == 0
    assert percent_A == 0

    seq = 'GCGCAGCAGGGGTGGTCGCCATGGAGACGCGTGGCCCTGGCCTGGCGGTCCGCGCTGAGAGTCGCCGATTAGTCGGCATCGGGCCTCGGGCGCCCCCGGGGCGGGTTGGGTTGCAGCCCAGCGGGCGGCTGGACCGCCGCGGTGGGGCGGGGACAATGGGGTACAAGGACAACGACGGCGAGGAGGAGGAGCGGGAGGGCGGCGCCGCGGGCCCGCGGGGGTCTAGACTGCCCCCCATCACAGGCGGCGCCTCCGAGCTGGCCAAACGGAAGGTGAAGAAGAAAAAAAGGAAGAAGAAGACCAAGGGGTCTGGCAAGGGGGACGGGATCTTGCTCTGGGGAAAGCCAGCTGCAGTGTTGTGTGGGCTGCCTGAGGAAAGGTGCCCTCCGGAAAGAAATGATGTCTCCGCCCAACAGCATGGTGGACCTGAGTCCTGCTACCAGCCACATGAGTGTGCTTGGAAGCAGGTGGAGCCCAGTTGAGCCTTGAGGTGCCTCCCTCTGTGAGAGACCCCAAGCCGGAGACATCCGGCCAAGAGGCACCCAGATTCCTGCCCCACAGAAACTGATACAATAACGTTGTTTAAAGCCACTACATTTGGAGGTAATCTTGTTACACAGCAATAACTGACTAATACAGTCGGCCTCTTCCATCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    cigar = '324M10025N331M100S'
    flag = '0'
    pos = '410323'
    read = load_reads(seq, cigar, flag, pos)

    polyA_site, tail_len, percent_A = PolyaTailCounter.detect_polyA_tail(read)

    assert polyA_site == 421002
    assert tail_len == 100
    assert percent_A == 1

    seq = 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATGACGTTTCCATTTACTTTGGATTATATGTCATTATAAATATTAACAAATAAGACTTAAAAAGGACACCTTCGGGTAGGTCAGACCAAAATACAAAACTTGTCTGTGGGACTGCAGTTTGAGGACAGTGTCTGCAGCCGTCACATGGCAGCAAAACGGGGTTAAGCAGTGCACGAGAGTCTGCGTCGACGACAGCCAGAGTCCATGCATCGGGAGGTTCACTCGGTTTGCGAAGGAACAACGGGCTCGGCATGCACGGGCTCGGCGGCGGCGGACGGGCCGGGGCGCAGTTCCCCGCGCTCGCCACTAGAGGTCAGGAGGTGACCGCTTCGGGGCTGGAAGACGGGCCCGTCGGGGATTGGCGCAGGCGGCGGGCGGGGCGGCGGGCGGGGCGGCGCTGGAGGCAGCGCCTGGTTACTGACACCTGGAATGACTTTTTTTTTTTGGCATCAGATTTCCTGTCTTTGTGGGGATGATGGACCCGAGTAAAGATGCCCATTCGGGGTCAAAGGCAGAGCCGCTTCTGCAGCTTCTCAAAGCGTTGTTTGTTTGTTTTTTTTTCTGAGACGGAGTCTTGCTCTGTCGCCCAGGCTGGAGTGCAGTGCCGCGATCTTGGCTCACTGCAGCCTTCACCTCCCGGGTTCAAGCGATTCTCGTGCCTCAGCCTCCCTGAGTAGCTGGGACTACAGGCGTGCGCCACCACACCCGGCTAATTTTTGTATTTTTAGTAGAGATGGGGTTTCGCCCTGTTGGCCAGGCGGGTTTCGAACTCCTGACCTCAGGTGATCCGCCAGCCTCGGCCTCCCAAAGTGCTGCAATTACAGGCGTGAGCCATCGTGCCTCAAAACGCTTTAACAGAAAGACAATCTGCACGGGATCTAAAAGGGTGCTGAGATCCTAGGGAAGGAAGGATCCAAACTTCCTGGGGAGTTCCTGCCCGAGTGCCTGTGCTGCCCCTGGGCTGGCTGGCCAGTAAGCCCGCCTCCCAGCCTGACTGTCCCCATCTTTCGGTCCCAGCCCCATTTGCACAGCCTGGGCAGTAGAGGGCCCTGGACTGGGGGGCTGGAGTTCTGGGTTCTTGTCCCAGCTGTGCCCCTTCAGGCCCCTTTCACCCACTGGGCCTCACATTCCCCATGTGCCTAATAAGGAAGTCATGGTAGGTGGGGGGTACAGCCTCTTCTCGCTTTGCCATTCACTTCTGTGTCTTCAGAGAGGCTGTCAAATCTCCTCTCTGAATGAGACCCCCAAAAAAGCAAAGCTAAGAAGATACCCAGAGCTTAACTAGACACCAGGCCTTTTAGAAATAGACCACCTCTTACCTTAGGCCCCCAGAGGGTGCCCATTCTGTGTGGAGAAAAGAGGAGCCCTTGCCTCAGCCCCCAGAGGCTAGGGTGGGGTGGCTGAGTTTTGGGGCCAGGTTAGACACCTCTGGGGAAGCCTGAAGTAGCACGGATGGTTTCAAAGCCAGCGGATGAGGTGGCAGGACAGTGACAAGACCCCAGGTCTCCATCATGCACCTGAGCTTACTGAGCCTTCACCTGGTGTCTCTGCTGAGTCTCCGACAGACCAAGAGGGAAGGGACTGGGGTACCGACCCCCAGAGAGAGAAAGCGGCAGTCAGAGGACCTGGGTTTGAGTCCTGGTTCACCCCTTCCTGGCTGTGTGGCCTTGGCAAACTACTCAGACTCTGGGAACCTGTTTCACCTGCAAGATGGGGATGAGAATCAAACCCACCTTGCAGGGCTGTGAACATGCGTTCAGACCAGTGCATGCAAAAACCTTCACACAAAAACCCTTCCTAAAGGCAGGGAGACAAGTCTCCAGGAGCAGCAGGTGGCCAGGGACTGTGTGGGGGCTGGGGTCCTGTTTTCCCCGCAACCTGGGAAAGGCCTGATGGGCACTTGGTAGGATTCAAATCATAACCCTGGACCTCAGGTGGTGGTGTGTGCTTTGTGTCTGCAGTGGAAGGTCCCACCGTCGTGCCTGTGAGTCCCTCTATCGTGTTGAAGGGGTGGGCAGCAGGCGGGGGAGCCCAGGCTGGCAGCTAGCAGGACCTCTCTGGTGTGAGCTCAGCACGCCAATTCCCCTGAAGCGTGGTGTCCAGCACTCTGGGCTGGGGGCTGTGGATCCTGGTTCCAGCTGTGTGGGGCCTGGAAGGCCCTGGGCAGGTCACCTGACCTCTCTGGGCCTGTTTCCATCACACACTGATGGGCTGAGCACACTGGGACTCTGGCTGTGCAACTCCTTGGCTCTGCATTTGTTCACCCAGCGTTCCTGAGGGGCCCTTGGTAGGCAGAGAAAGTTTGTGGGCTTCAGGCATGGGCCCCAAGATTCAGATACTCTCAAGCCTCCTGGGGAGTCTCACTCAGGGGAGCAGACAGGCCCACCAGCCAGGGTGATTCTTGCTCAGCTCCTGGAGGGTGGAGCTGGCCCCACAGGCCTCCCCAAAGACAAGGCCCTGGACGGCCACTGATTACTCCCAAAAGGGACATCTGTGGCGGTAGTGGGGACCAAGATGCCACGAGCAAGCACCCGGAGGCCCCCAGTGTGAGCCATGGAGTGGAGGGGAGGGGAAAGGGCAGAGTCAGGACGTGTAGGAATGCTTGCTTTTTTTCCAAGCACAAGGGACCCTTTTCTCCACTGCAGCTGACCTGATGCTTATGCCAGGAAGGAGGAGGGGCGGGCTCCGTCCTTGAGGTCCCTCAGGAGTAGAAAGAGATCAGAGTGGGAGACTTGGGTCTGAGGTCTGAACTTGAGCCTACACCAGTTTCTCCATGGTGTTTCCATCACGTCCCCCACCTCCTGCCTCGAGCCTCACACCTTCCTAACAAACCCTCCCCTGGAGAGGAGACCCTGGGTCCACAGCACCCGGCCCCATTGGCTTTCTCTCTCCAGGCGTTTTAGACCACTCCCACTGCCTGGACTGCTTTTAGAAGCCCTCACTCACCCTCCAAGGCCTTACGGAAGCACCGCCTCCTCCAGGAAGCCGTCCCTGACCTCCCGGCAGAGTCAGCAGAGCCCACGGTACTTGCAGACTGGCCAAGCTGGTCTTGGTATTTCTCACCCCAGTTGGAATGGTTTGTTCCCAGCTTCCTGACTAGATGGGAACTCCCTGAGGGCAGGCCCTGTGTCTCATTCACCCCAGGGCTTGTGGAATCATTGAGTGAAGCATGTGCCAATTTACCCATCATCAGGAGCCTCCGGGAACTCTGGCAGACTTTCGGCCGGCAGGCCCGCTTCTTCCATCTGTCCAGCAGCGAAGGAGATGAGGGATGCAGTTAGGCTTTCTTGGGCTGGAGCAGCCAGTCTTCAAGGTCCCATCCTCCACC'
    cigar = '104S257M5D11M2I97M12D186M1I1030M1D357M2I57M1D31M1D474M1I108M1D253M1D272M1D190M'
    flag = '16'
    pos = '143699'
    read = load_reads(seq, cigar, flag, pos)

    polyA_site, tail_len, percent_A = PolyaTailCounter.detect_polyA_tail(read)

    assert polyA_site == 143698
    assert tail_len == 104
    assert percent_A == 1

    seq = 'CAAATAGGAGCAGTATGATTTTTTCTTTTATTTTTAATTGACAAATAATAATTATATAAACAATGTAATTTCTAACAACAAAGGAATGCTAAAAAATATACAAAAAAAAAAAAAAAAAAAAGAGCGGAGGCGCCCCCGGGTGCACGCCGG'
    cigar = '1S119M30S'
    flag = '0'
    pos = '716075'
    read = load_reads(seq, cigar, flag, pos)

    polyA_site, tail_len, percent_A = PolyaTailCounter.detect_polyA_tail(read)

    assert polyA_site == 716193
    assert tail_len == 1
    assert percent_A == 1

    seq = 'GGGCCGGGCGGGCCGCCGGGGGCGCTGCCGCCCTTTTTTTTTTTTTTTTTTAAGAGCAGAGTGAACTCTTTATTGATTATACAAATTACCACTATTTATTTTAAACCCAAAGTGACTTCAAAGTTGTTTTTTGGTTTTTAAAGGGGCTAC'
    flag = '16'
    cigar = '50S100M'
    pos = '515042'
    read = load_reads(seq, cigar, flag, pos)

    polyA_site, tail_len, percent_A = PolyaTailCounter.detect_polyA_tail(read)

    assert polyA_site == 515041
    assert tail_len == 17
    assert percent_A == 1


@pytest.fixture
def tail_counter_short():
    return PolyaTailCounter(short_bam)


def test_PolyaTailCounter_iter_tailed_reads(tail_counter_short):
    tailed = sum(1 for i in tail_counter_short.iter_tailed_reads())
    assert tailed > 3900


def test_PolyaTailCounter_tail_len_dist(tail_counter_short):
    tail_len = tail_counter_short.tail_len_dist()
    assert tail_len.shape[0] > 20


def test_PolyaTailCounter_save_tailed_reads(tail_counter_short, tmp_path):
    tailed_bam = str(tmp_path / 'tailed.bam')
    tail_counter_short.save_tailed_reads(tailed_bam)

    tailed_bam = pysam.AlignmentFile(tailed_bam, 'rb')
    tailed = sum(1 for i in tailed_bam)
    assert tailed > 3900


def test_PolyaTailCounter_count(tail_counter_short):
    tes = tail_counter_short.count()
    assert len(tes) > 1000


def test_PolyaTailCounter_to_df(tail_counter_short):
    df = tail_counter_short.to_df()
    assert df.shape[0] > 1000


def test_TesMultiCounter_samples(tmp_path):
    # Test multiple samples
    df = pd.DataFrame({
        'sample': ['short_rep1', 'short_rep2', 'short', 'short'],
        'path': [short_bam, short_bam, short_bam, short_bam]
    })
    csv_path = tmp_path / 'samples.csv'
    df.to_csv(csv_path)

    counter = TesMultiCounter(csv_path, method='tail')
    df_all, tes = counter.to_df()
   
    pd.testing.assert_frame_equal(
        tes['short_rep1'],
        tes['short_rep2']
    )
    _df = tes['short_rep1'].copy()
    _df['count'] *= 4
    pd.testing.assert_frame_equal(df_all, _df)

    _df = tes['short_rep1'].copy()
    _df['count'] *= 2
    pd.testing.assert_frame_equal(tes['short'], _df)

    # Test bigwig exporting
    output_dir = tmp_path / 'multi_sample'
    output_dir.mkdir()    
    TesMultiCounter._to_bigwig(df_all, tes, chr17_chrom_sizes,
                               output_dir, prefix='tes_counts')

    assert {i.name for i in output_dir.iterdir()} == {
        'all_tes_counts_pos.bw',
        'all_tes_counts_neg.bw',
        'short_tes_counts_pos.bw',
        'short_tes_counts_neg.bw',
        'short_rep1_tes_counts_pos.bw',
        'short_rep1_tes_counts_neg.bw',
        'short_rep2_tes_counts_pos.bw',
        'short_rep2_tes_counts_neg.bw'
    }

    bw_pos = pyBigWig.open(str(output_dir / 'all_tes_counts_pos.bw'))
    bw_neg = pyBigWig.open(str(output_dir / 'all_tes_counts_neg.bw'))
    count = np.nansum(
        bw_pos.values('chr17',
                      int(df_all['End'].min() - 1),
                      int(df_all['End'].max()))
    ) + np.nansum(
        bw_neg.values('chr17',
                      int(df_all['End'].min() - 1),
                      int(df_all['End'].max()))
    )
    assert df_all['count'].sum() == count

    bw_pos = pyBigWig.open(str(output_dir / 'short_tes_counts_pos.bw'))
    bw_neg = pyBigWig.open(str(output_dir / 'short_tes_counts_neg.bw'))
    count = np.nansum(
        bw_pos.values('chr17',
                      int(df_all['End'].min() - 1),
                      int(df_all['End'].max()))
    ) + np.nansum(
        bw_neg.values('chr17',
                      int(df_all['End'].min() - 1),
                      int(df_all['End'].max()))
    )
    assert df_all['count'].sum() / 2 == count

    # Test single sample
    df = pd.DataFrame({
        'sample': ['short', 'short'],
        'path': [short_bam, short_bam]
    })
    df.to_csv(csv_path)

    counter = TesMultiCounter(csv_path, method='tail')
    df_all, tes = counter.to_df()

    output_dir = tmp_path / 'single_sample'
    output_dir.mkdir()
    TesMultiCounter._to_bigwig(df_all, tes, chr17_chrom_sizes,
                               output_dir, prefix='tes_counts')
    assert {i.name for i in output_dir.iterdir()} == {
        'all_tes_counts_pos.bw', 'all_tes_counts_neg.bw'}

    # Test read_annot
    counter = TesMultiCounter(read_annot)
    df_all, tes = counter.to_df()

    output_dir = tmp_path / 'read_annot_dir'
    output_dir.mkdir()
    TesMultiCounter._to_bigwig(df_all, tes, chr17_chrom_sizes,
                               output_dir, prefix='tes_counts')

    assert {i.name for i in output_dir.iterdir()} == {
        'all_tes_counts_pos.bw',
        'all_tes_counts_neg.bw',        
        'gm12878_tes_counts_neg.bw',
        'gm12878_tes_counts_pos.bw',
        'hepg2_tes_counts_neg.bw',
        'hepg2_tes_counts_pos.bw'
    }
    


def test_TssMultiCounter_samples(tmp_path):
    # Test read_annot
    counter = TssMultiCounter(read_annot)
    df_all, tes = counter.to_df()

    output_dir = tmp_path / 'read_annot_dir'
    output_dir.mkdir()
    TssMultiCounter._to_bigwig(df_all, tes, chr17_chrom_sizes,
                               output_dir, prefix='tss_counts')

    assert {i.name for i in output_dir.iterdir()} == {
        'all_tss_counts_pos.bw',
        'all_tss_counts_neg.bw',
        'gm12878_tss_counts_neg.bw',
        'gm12878_tss_counts_pos.bw',
        'hepg2_tss_counts_neg.bw',
        'hepg2_tss_counts_pos.bw'
    }
