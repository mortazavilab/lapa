bam = snakemake.params['bam']
sample = snakemake.wildcards['sample']

with open(snakemake.output['csv'], 'w') as f:
    f.write('sample,path\n')

    for i in ['1', '2', '3']:
        f.write(f'{sample}_{i},{bam.format(sample=sample, rep=i)}\n')
