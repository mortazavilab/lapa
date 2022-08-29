import pysam


bam = pysam.AlignmentFile(snakemake.input['bam'])
stats = open(snakemake.output['txt'], 'w')

for read in bam:
    if read.mapping_quality >= 10 and (not read.is_secondary):
        segments = read.get_aligned_pairs(matches_only=True)
        stats.write('%d\n' % len(segments))
