# from itertools import chain
# import numpy as np
# import pandas as pd
# import pyranges as pr
# from tqdm import tqdm
# from lapa.utils.io import read_polyA_cluster
# from lapa.cluster import TesClustering


# def talon_tes_read_mapping(df_read_annot, fasta):
#     df_read_annot['End'] = np.where(df_read_annot['Strand'] == '-',
#                                     df_read_annot['Start'],
#                                     df_read_annot['End'])
#     del df_read_annot['Start']
#     df_read_annot['count'] = 1

#     df_read_annot = df_read_annot.rename(columns={'read_name': 'reads'})
#     df_read_annot = df_read_annot.groupby(
#         ['Chromosome', 'End', 'Strand', 'gene_id']).agg(
#         {'reads': list, 'count': 'sum'}).reset_index()
#     clustering = TesClustering(fasta, groupby='gene_id', fields=[
#                                'reads'], progress=False)

#     df = list()

#     for _, _df in tqdm(df_read_annot.groupby('gene_id')):
#         clusters = list(clustering.cluster(_df))
#         polyA_site = {c: c.polyA_site() for c in clusters}
#         # TO FIX: filter minimum count, internal priming

#         mapping = dict()
#         for c in clusters:
#             for reads in c.fields['reads']:
#                 for r in reads:
#                     mapping[r] = polyA_site[c]

#         for _, row in _df.iterrows():
#             for read in row['reads']:
#                 if read in mapping:
#                     df.append((read, mapping[read]))
#                 else:
#                     df.append((read, -1))

#     return pd.DataFrame(df, columns=['read_name', 'polyA_site'])
