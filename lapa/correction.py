import pandas as pd
import pyranges as pr
from tqdm import tqdm
from lapa.utils.io import read_talon_read_annot


class Transcript:
    '''
    Transcript class for performing manupulation on transcripts.

    Args:
      transcript_id: transcript id
      df: Transcript and subfeatures as data.frame. 
        Enteries of gtf file related with the transcript.
      min_exon_len: Minimum exon length
    '''

    def __init__(self, transcript_id, df, min_exon_len=25):
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
        '''
        Most five prime exon of transcript
        '''
        if self.strand == '+':
            return self.first_exon_idx
        elif self.strand == '-':
            return self.last_exon_idx

    @property
    def three_prime_exon_idx(self):
        '''
        Most three prime exon of transcript
        '''        
        if self.strand == '+':
            return self.last_exon_idx
        elif self.strand == '-':
            return self.first_exon_idx

    def copy(self, new_transcript_id):
        '''
        Create copy of transcript with new transcript id.

        Args:
          new_transcript_id: `oldTranscriptId#suffix`
        '''
        old_transcript, prefix = new_transcript_id.split('#')
        assert old_transcript == self.transcript_id

        df = self.df.copy()
        df['transcript_id'] = new_transcript_id
        df['exon_id'] = df['exon_id'] + '#' + prefix
        return Transcript(new_transcript_id, df)

    def valid_five_prime_exon_len(self, tss_site: int):
        '''
        Checks new proposed tss_site is valid for most five
        prime exon of transcript based on coordinates and 
        minimum exon length.

        Args:
          tss_site: Position of proposed tss site 
        '''
        exon_idx = self.five_prime_exon_idx

        if self.strand == '+':
            exon_len = self.df.loc[exon_idx, 'End'] - tss_site
        elif self.strand == '-':
            exon_len = tss_site - self.df.loc[exon_idx, 'Start']

        return exon_len > self.min_exon_len

    def valid_three_prime_exon_len(self, polyA_site):
        '''
        Checks new proposed polyA_site is valid for most three 
        prime exon of transcript based on coordinates and 
        minimum exon length.

        Args:
          polyA_site: Position of proposed poly(A) site
        '''
        exon_idx = self.three_prime_exon_idx

        if self.strand == '+':
            exon_len = polyA_site - self.df.loc[exon_idx, 'Start']
        elif self.strand == '-':
            exon_len = self.df.loc[exon_idx, 'End'] - polyA_site

        return exon_len > self.min_exon_len

    def update_tss_site(self, tss_site):
        '''
        Updates tss site of transcript and most five prime exons.

        Args:
          tss_site: Position of proposed tss site.        
        '''
        exon_idx = self.five_prime_exon_idx

        if not self.valid_five_prime_exon_len(tss_site):
            raise ValueError(
                f'Exon length is shorter than min_exon_len={self.min_exon_len}'
                f' for transcript={self.transcript_id} '
                f' and tss_site={tss_site}')

        if self.strand == '+':
            self.df.loc[exon_idx, 'Start'] = tss_site
            self.df['Start'].iloc[0] = tss_site

        if self.strand == '-':
            self.df.loc[exon_idx, 'End'] = tss_site
            self.df['End'].iloc[0] = tss_site

    def update_polyA_site(self, polyA_site):
        '''
        Updates poly(A) site of transcript and most three prime exons.

        Args:
          polyA_site: Position of proposed poly(A) site.
        '''        
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
    '''
    Modifier to update transcript start, end sites
    of transcript and respective exons, genes.

    Args:
      templete_gtf: Use gtf as templete.
      min_exon_len: Minimum exon length.
    '''

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
        '''
        Fetch transcript from the templete gtf and return transcript object.
        '''
        return Transcript(
            transcript_id,
            self._transcript_templetes[transcript_id],
            min_exon_len=self.min_exon_len
        )

    def add_transcript(self, transcript):
        '''
        Add new trascript isoform to modifier.
        '''        
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
        '''
        Save all motifiers with motified trascript as gtf.

        Args:
          path: Output path to save gtf.
        '''
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
        'transcript_id', 'Strand', 'tss_site', 'polyA_site', 'sample'
    ]).agg({
        'read_Start': 'min', 'read_End': 'max', 'count': 'sum'
    }).reset_index()


def _transcript_tss_tes(df, threshold=1):
    '''
    '''
    df = df[(df['polyA_site'] != -1) & (df['tss_site'] != -1)]
    df = df.set_index(['transcript_id', 'tss_site', 'polyA_site'])

    _df = df.groupby(['transcript_id', 'tss_site',
                      'polyA_site']).agg({'count': 'sum'})
    df = df[_df['count'] > threshold]

    # update transcript_ids
    _df = df.reset_index().drop_duplicates(
        subset=['transcript_id', 'tss_site', 'polyA_site'])
    _df['suffix'] = _df.groupby('transcript_id') \
                       .cumcount().astype(str).radd('#')
    df = df.join(_df.set_index(
        ['transcript_id', 'tss_site', 'polyA_site'])['suffix']).reset_index()
    df['transcript_id'] += df['suffix']
    del df['suffix']

    return df[['transcript_id', 'tss_site', 'polyA_site', 'sample', 'count']]


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

        if not transcript.valid_five_prime_exon_len(int(row.tss_site)):
            print(f'Transcript {row.transcript_id} is not corrected '
                  'because exon chain longer than ends')
            continue
        transcript.update_tss_site(int(row.tss_site))

        if not transcript.valid_three_prime_exon_len(int(row.polyA_site)):
            print(f'Transcript {row.transcript_id} is not corrected '
                  'because exon chain longer than ends')
            continue
        transcript.update_polyA_site(int(row.polyA_site))

        modifier.add_transcript(transcript)

    modifier.to_gtf(gtf_output)


def _update_abundace(df_abundance, df_link_counts, keep_unsupported=False):

    cols_abundance = df_abundance.columns[:11]
    samples_abundance = df_abundance.columns[11:]

    if set(df_link_counts['sample'].unique()) != set(samples_abundance):
        raise ValueError(
            'Samples in abundance file do not match with read_annot file. '
            f'read_annot samples: {set(df_link_counts["sample"].unique())} '
            f'abundance samples: {set(samples_abundance)} ')

    df = df_link_counts.pivot(index='transcript_id',
                              columns='sample', values='count') \
        .fillna(0).astype(int).reset_index()

    df['annot_transcript_id'] = df['transcript_id'].str.split('#').str.get(0)

    df_abundance_lapa = df \
        .set_index('annot_transcript_id') \
        .join(df_abundance[cols_abundance].set_index('annot_transcript_id'),
              how='inner') \
        .reset_index(drop=True) \
        .rename(columns={'transcript_id': 'annot_transcript_id'}) \
        .set_index('annot_transcript_id')

    df_abundance_lapa['annot_transcript_name'] += '#' \
        + df_abundance_lapa.index.str.split('#').str.get(1)

    if keep_unsupported:
        df = df_abundance_lapa.reset_index()
        df['annot_transcript_id'] = df['annot_transcript_id'] \
            .str.split('#').str.get(0)
        df = df.groupby('annot_transcript_id').sum()

        df_abundance = df_abundance.set_index('annot_transcript_id')
        transcripts = list(set(df_abundance.index).intersection(df.index))

        for i in samples_abundance:
            df_abundance.loc[transcripts, i] -= df.loc[transcripts, i]

        df_abundance_lapa = pd.concat([
            df_abundance,
            df_abundance_lapa
        ])

    return df_abundance_lapa.reset_index() \
        .sort_values('annot_transcript_id')


def correct_talon(links_path, read_annot_path, gtf_input,
                  gtf_output, abundance_path, abundance_output,
                  link_threshold=1, keep_unsupported=False):
    '''
    LAPA creates GTF file with tss/poly(A) cluster support based
    on the linking reads and using splice chain of TALON.

    Args:
      links_path: Path to linking read file generated 
        with `lapa_link_tss_to_tes` command.

      read_annot_path: read_annot of TALON annotating read
        transcript assignments.
      gtf_input: Input gtf file to extract splice chains.
      gtf_output: Output corrected gtf contains trascripts
        with tss/poly(A) end support.
      abundance_path: Input abundance file of TALON which contains
        abundance of each transcript.
      abundance_output: Update abundance file which calculated
        based on abundance of linking reads.
      link_threshold: Minimum number of linking reads to create 
        transcript isoform.
      keep_unsupported: Keep transcripts without tss and tes support in
        the original gtf. If true transcript created with
        non-linking reads (partial) in the original files
        are kept gtf and abundance.
    '''
    df = _links_transcript_agg(links_path, read_annot_path)
    df = _transcript_tss_tes(df, threshold=link_threshold)

    df_abundance = pd.read_csv(abundance_path, sep='\t')
    df_abundance_cor = _update_abundace(df_abundance, df, keep_unsupported)
    df_abundance_cor[df_abundance.columns].to_csv(
        abundance_output, index=False, sep='\t')

    df = df.drop_duplicates(
        subset=['transcript_id', 'tss_site', 'polyA_site'])
    _save_corrected_gtf(df, gtf_input, gtf_output, keep_unsupported)
