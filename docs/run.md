## Poly(A) cluster calling

Poly(A) clusters can be created with followng line:

```
lapa --alignment {bam} --fasta {fasta} --annotation {gtf} --chrom_sizes {chrom_sizes} --output_dir {output}
```

Output of LAPA (`--output_dir`) is directory and content explained in output section of documentation.

Fasta is standart genome sequence file can be downloaded (Ensembl or GENCODE). Based on the FASTA file `chrom_sizes` file can be generated with the following command:

```
faidx {fasta} -i chromsizes > {chrom_sizes}
```

GTF file similarly can be downloaded from (Ensembl or GENCODE). Ensembl separately annotates UTR features of GTF file as **three_prime_utr** and **five_prime_utr** so directly used as input to LAPA. Yet GENCODE does not distinguishe 5' UTR or 3' UTR and annotated all as just **UTR**; thus, UTR entry of GENCODE need to be corrected before passing it to LAPA with following command with tool called `gencode_utr_fix`:

```
pip install cython
pip install -e git+https://github.com/MuhammedHasan/gencode_utr_fix.git#egg=gencode_utr_fix

gencode_utr_fix --input_gtf {gencode.gtf} --output_gtf {gencode.utr_fixed.gtf}
```

LAPA can process single or multiple bam files as input. For example, single bam file:

```
lapa --alignment single.bam --fasta {fasta} --annotation {gtf} --chrom_sizes {chrom_sizes} --output_dir {output}
```

> **_Attention!:_** Running LAPA with multiple samples joinly is high prefered over running it on individuals files separately.

Multiple bam files can be given as input with comma separation:

```
lapa --alignment sample1.bam,sample2.bam,sample3.bam  --fasta {fasta} --annotation {gtf} --chrom_sizes {chrom_sizes} --output_dir {output}
```

If multiple bam files are given as input LAPA filters poly(A) clusters based on the replication rate between samples. (See the documentation for CLI). If comma separated bam files provided with the command above, LAPA assumes all the files are experimental replicates of each other.

For multiple samples, csv file can also be provided as alignment input then replication rate will be calculated based on the dataset column:

`samples.csv`:

```
sample,dataset,path
biosample1_rep1,tissue1,bam_sample1_rep1_path.bam
biosample1_rep2,tissue1,bam_sample1_rep2_path.bam
biosample2_rep1,tissue2,bam_sample2_rep1_path.bam
biosample2_rep2,tissue2,bam_sample2_rep2_path.bam
```

Then:

```
lapa --alignment samples.csv  --fasta {fasta} --annotation {gtf} --chrom_sizes {chrom_sizes} --output_dir {output}
```


Replication rate is (in my opinion as the developer of LAPA) most important parameter of this package because it determines the final set of poly(A) clusters. If you have experimental replicates pass them as dataset in the sample file. By default, LAPA assumes datasets are experimental replicates of each other and applies conservative 0.95 replication rate cutoff. If experimental replicates rather biological replicates with batch effect .etc default replication rate could be stringent. For higher recall setting, replication rate  could be set to lower cutoff such as (`--min_replication_rate 0.75`). Moreover, default replication rate accepts clusters as replicated if there is at least 1 read 2 out of N samples. If you have many ~10-100 replicates, you may set this to `replication_num_sample` to higher number just of 5 or 10% of samples. To visualize and choice right replication rate, see the quality control section.


## TSS cluster calling

Options for TSS cluster calling for similar so see [the documentation of Poly(A) clusters](#poly(A)-clusters-calling) section.

```
lapa_tss --alignment sample1.bam,sample2.bam,sample3.bam  --fasta {fasta} --annotation {gtf} --chrom_sizes {chrom_sizes} --output_dir {output}
```

Or:


```
lapa --alignment samples.csv  --fasta {fasta} --annotation {gtf} --chrom_sizes {chrom_sizes} --output_dir {output}
```


## Linking reads

After performing TSS and poly(A) cluster calling, linking reads can be calculated from long-reads. Linking reads are the reads start in TSS cluster and ends in poly(A) cluster. 

```
lapa_link_tss_to_tes --alignment sample_rep1.bam --lapa_dir {lapa_dir} --lapa_tss_dir {lapa_tss_dir} --output {linking_reads}.csv
```

If dataset specific linking read calls can be obtained with:


```
lapa_link_tss_to_tes --alignment sample_rep1.bam --lapa_dir {lapa_dir} --lapa_tss_dir {lapa_tss_dir} --output {linking_reads}.csv --dataset {dataset_name}
```


## Correction of TALON gtf and abundance


```
lapa_correct_talon \
    --links {links_reads}.csv \
    --read_annot {talon_read_annot} \
    --gtf_input {talon_gtf} \
    --gtf_output {corrected_gtf} \
    --abundance_input {talon_abundance} \
    --abundance_output {corrected_abundance}
```
