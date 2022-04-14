import pandas as pd
import pyranges as pr
from tqdm import tqdm
from lapa.utils.io import read_talon_read_annot


class Transcript:

    def __init__(self, transcript_id, df, min_exon_len=25):
        '''
        Args:
          df: enteries of gtf file related with the transcript
        '''
        assert df['Feature'].iloc[0] == 'transcript', \
            "The first features of df should be transcript"

        self.transcript_id = transcript_id
        self.df = df
        self.min_exon_len = min_exon_len

        self.gene_id = df['gene_id'].iloc[0]
        self.strand = df['Strand'].iloc[0]
        assert self.strand in {'+', '-'}, \
            f"Strand should be one of `+` or `-` but got {self.strand}"

        _df_exon = df[df['Feature'] == 'exon']
        self.first_exon_idx = _df_exon['Start'].idxmin()
        self.last_exon_idx = _df_exon['End'].idxmax()

    @property
    def five_prime_exon_idx(self):
        if self.strand == '+':
            return self.first_exon_idx
        elif self.strand == '-':
            return self.last_exon_idx

    @property
    def three_prime_exon_idx(self):
        if self.strand == '+':
            return self.last_exon_idx
        elif self.strand == '-':
            return self.first_exon_idx

    def copy(self, new_transcript_id):
        '''
        Args:
          new_transcript_id: `oldTranscriptId_suffix`
        '''
        old_transcript, prefix = new_transcript_id.split('#')
        assert old_transcript == self.transcript_id

        df = self.df.copy()
        df['transcript_id'] = new_transcript_id
        df['exon_id'] = df['exon_id'] + '#' + prefix
        return Transcript(new_transcript_id, df)

    def valid_five_prime_exon_len(self, start_site):
        exon_idx = self.five_prime_exon_idx

        if self.strand == '+':
            exon_len = self.df.loc[exon_idx, 'End'] - start_site
        elif self.strand == '-':
            exon_len = start_site - self.df.loc[exon_idx, 'Start']

        return exon_len > self.min_exon_len

    def valid_three_prime_exon_len(self, polyA_site):
        exon_idx = self.three_prime_exon_idx

        if self.strand == '+':
            exon_len = polyA_site - self.df.loc[exon_idx, 'Start']
        elif self.strand == '-':
            exon_len = self.df.loc[exon_idx, 'End'] - polyA_site

        return exon_len > self.min_exon_len

    def update_start_site(self, start_site):
        exon_idx = self.five_prime_exon_idx

        if not self.valid_five_prime_exon_len(start_site):
            raise ValueError(
                f'Exon length is shorter than min_exon_len={self.min_exon_len}'
                f' for transcript={self.transcript_id} '
                f' and start_site={start_site}')

        if self.strand == '+':
            self.df.loc[exon_idx, 'Start'] = start_site
            self.df['Start'].iloc[0] = start_site

        if self.strand == '-':
            self.df.loc[exon_idx, 'End'] = start_site
            self.df['End'].iloc[0] = start_site

    def update_polyA_site(self, polyA_site):
        exon_idx = self.three_prime_exon_idx

        if not self.valid_three_prime_exon_len(polyA_site):
            raise ValueError(
                f'Exon length is shorter than min_exon_len={self.min_exon_len}'
                f' for transcript={self.transcript_id}'
                f' and polyA_site={polyA_site}')

        if self.strand == '+':
            self.df.loc[exon_idx, 'End'] = polyA_site
            self.df['End'].iloc[0] = polyA_site

        if self.strand == '-':
            self.df.loc[exon_idx, 'Start'] = polyA_site
            self.df['Start'].iloc[0] = polyA_site


class TranscriptModifier:

    def __init__(self, templete_gtf, min_exon_len=25):
        self.gtf = pr.read_gtf(templete_gtf).df
        self.min_exon_len = min_exon_len

        _groups = self.gtf[~self.gtf.transcript_id.isna()].groupby(
            'transcript_id')
        self._transcript_templetes = {
            transcript_id: df
            for transcript_id, df in _groups
        }

        _groups = self.gtf[(self.gtf['Feature'] == 'gene')].groupby('gene_id')
        self._gene_templete = {
            gene_id: df
            for gene_id, df in _groups
        }

        self.transcripts = dict()
        self.genes = dict()

    def fetch_transcript(self, transcript_id):
        return Transcript(
            transcript_id,
            self._transcript_templetes[transcript_id],
            min_exon_len=self.min_exon_len
        )

    def add_transcript(self, transcript):

        self.transcripts[transcript.transcript_id] = transcript.df

        if transcript.gene_id not in self.genes:
            self.genes[transcript.gene_id] = self._gene_templete[
                transcript.gene_id]

        df = self.genes[transcript.gene_id]

        df['Start'].iloc[0] = min(
            df['Start'].iloc[0],
            transcript.df['Start'].iloc[0])

        df['End'].iloc[0] = max(
            df['End'].iloc[0],
            transcript.df['End'].iloc[0])

        self.genes[transcript.gene_id] = df

    @staticmethod
    def _sort_gtf_key(col):
        if col.name == 'End':
            return -col
        elif col.name == 'exon_number':
            return col.astype(float)
        else:
            return col

    @staticmethod
    def _sort_gtf(df):
        return df.sort_values([
            'Chromosome', 'gene_id', 'transcript_id', 'exon_number',
            'Start', 'End'
        ], na_position='first', key=TranscriptModifier._sort_gtf_key)

    def to_gtf(self, path):
        print('Sorting and writing gtf...')
        df_gtf = pd.concat([
            *self.genes.values(),
            *self.transcripts.values()
        ])
        df_gtf = TranscriptModifier._sort_gtf(df_gtf)
        pr.PyRanges(df_gtf).to_gtf(path)

    def __contains__(self, transcript_id):
        return transcript_id in self._transcript_templetes


def _links_transcript_agg(links, read_annot_path):
    df_links = pd.read_csv(links).set_index('read_name')
    df_read_annot = read_talon_read_annot(read_annot_path)[
        ['sample', 'read_name', 'transcript_id']].set_index('read_name')

    df = df_links.join(df_read_annot, how='inner')

    df['count'] = 1
    return df.groupby([
        'transcript_id', 'Strand', 'start_site', 'polyA_site', 'sample'
    ]).agg({
        'read_Start': 'min', 'read_End': 'max', 'count': 'sum'
    }).reset_index()


def _transcript_tss_tes(df, threshold=1):
    '''
    '''
    df = df[(df['polyA_site'] != -1) & (df['start_site'] != -1)]
    df = df.set_index(['transcript_id', 'start_site', 'polyA_site'])

    _df = df.groupby(['transcript_id', 'start_site',
                     'polyA_site']).agg({'count': 'sum'})
    df = df[_df['count'] > threshold]

    # update transcript_ids
    _df = df.reset_index().drop_duplicates(
        subset=['transcript_id', 'start_site', 'polyA_site'])
    _df['suffix'] = _df.groupby(
        'transcript_id').cumcount().astype(str).radd('#')
    df = df.join(_df.set_index(
        ['transcript_id', 'start_site', 'polyA_site'])['suffix']).reset_index()
    df['transcript_id'] += df['suffix']
    del df['suffix']

    return df[['transcript_id', 'start_site', 'polyA_site', 'sample', 'count']]


def _save_corrected_gtf(df, gtf, gtf_output, keep_unsupported=False):
    '''
    '''
    modifier = TranscriptModifier(gtf)

    if keep_unsupported:
        for transcript_id in tqdm(modifier._transcript_templetes):
            transcript = modifier \
                .fetch_transcript(transcript_id)

            modifier.add_transcript(transcript)

    for row in tqdm(df.itertuples(), total=df.shape[0]):
        templete_transcript = row.transcript_id.split('#')[0]

        if templete_transcript not in modifier:
            continue

        transcript = modifier \
            .fetch_transcript(templete_transcript) \
            .copy(row.transcript_id)

        if not transcript.valid_five_prime_exon_len(int(row.start_site)):
            print(f'Transcript {row.transcript_id} is not corrected '
                  'because exon chain longer than ends')
            continue
        transcript.update_start_site(int(row.start_site))

        if not transcript.valid_three_prime_exon_len(int(row.polyA_site)):
            print(f'Transcript {row.transcript_id} is not corrected '
                  'because exon chain longer than ends')
            continue
        transcript.update_polyA_site(int(row.polyA_site))

        modifier.add_transcript(transcript)

    modifier.to_gtf(gtf_output)


def _update_abundace(df_abundance, df_link_counts):

    cols_abundance = df_abundance.columns[:12]
    samples_abundance = df_abundance.colums[12:]

    transcripts = set()

    # if samples_abundance check sample overlap

    df = df_link_counts.pivot(index='transcript_id',
                              columns='sample', values='count') \
        .fillna(0).astype(int).reset_index()
    df['annot_transcript_id'] = df['transcript_id'].str.split('#').str.get(0)

    __import__("pdb").set_trace()

    df_abundance_lapa = df.set_index('annot_transcript_id').join(
        df_abundance[cols_abundance].set_index('annot_transcript_id'),
        how='inner')

    df_abundance = df_abundance.set_index('annot_transcript_id')

    for i in samples_abundance:
        # __import__("pdb").set_trace()
        df_abundance[i] -= df.loc[df_abundance.index, i]

    return pd.concat([
        df_abundance,
        df_abundance_lapa
    ]).sort_values('annot_transcript_id')


def correct_talon(links_path, read_annot_path, gtf_input,
                  gtf_output, abundance_path=None,
                  link_threshold=1, keep_unsupported=False):
    '''
    '''
    df = _links_transcript_agg(links_path, read_annot_path)

    df = _transcript_tss_tes(df, threshold=link_threshold)

    df_abundance = pd.read_csv(abundance_path, sep='\t')

    _update_abundace(df_abundance, df)

    # subset

    _save_corrected_gtf(df, gtf_input, gtf_output, keep_unsupported)