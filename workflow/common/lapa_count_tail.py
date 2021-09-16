from lapa.count import TailTesCounter

count = sum(1 for _ in TailTesCounter(snakemake.input['bam']).iter_tailed_reads())

with open(snakemake.output['txt'], 'w') as f:
    f.write(str(count))
