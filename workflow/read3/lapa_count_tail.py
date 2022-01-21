from lapa.count import TailTesCounter

count = sum(
    1 for _ in TailTesCounter(
        snakemake.input['bam'],
        min_tail_len=2,
        min_percent_a=1
    ).iter_tailed_reads()
)

with open(snakemake.output['txt'], 'w') as f:
    f.write(str(count))
