# from conftest import read_annot, fasta, gtf, chrom_sizes
# from lapa.utils.io import read_talon_read_annot
# from lapa.talon import talon_tes_read_mapping


# def test_polyA_site_long_read(tmp_path):
#     df_read_annot = read_talon_read_annot(read_annot)
#     df = talon_tes_read_mapping(df_read_annot, fasta)

#     assert df_read_annot.shape[0] == df.shape[0]

#     assert df[df['polyA_site'] == -1].shape[0] < df_read_annot.shape[0] * 0.8

#     _df = df.set_index('read_name').join(df_read_annot.set_index('read_name'))

#     without_site = sum(_df.groupby('gene_id').agg(
#         {'polyA_site': 'max'})['polyA_site'] == -1)

#     import pdb
#     pdb.set_trace()
