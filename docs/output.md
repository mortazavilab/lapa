# Output format


## Poly(A) clusters

Output directory of LAPA (`lapa --output_dir your_output_dir`) looks like:

```
your_output_dir/
├── polyA_clusters.bed
├── raw_polyA_clusters.bed
├── counts
│   ├── all_polyA_counts_neg.bw
│   ├── all_polyA_counts_pos.bw
│   ├── {sample}_polyA_counts_neg.bw
│   ├── {sample}_polyA_counts_pos.bw
├── coverage
│   ├── all_polyA_coverage_neg.bw
│   ├── all_polyA_coverage_pos.bw
│   ├── {sample}_polyA_coverage_neg.bw
│   ├── {sample}_polyA_coverage_pos.bw
├── ratio
│   ├── all_polyA_ratio_neg.bw
│   ├── all_polyA_ratio_pos.bw
│   ├── {sample}_polyA_ratio_neg.bw
│   ├── {sample}_polyA_ratio_pos.bw
├── dataset
│   └── {dataset}.bed
├── raw_sample
│   └── {sample}.bed
├── sample
│   └── {sample}.bed
├── logs
│   ├── final_stats.log
└── ├── progress.log
    └── warnings.log
```

where `polyA_clusters.bed` is the main output of LAPA and contains poly(A) clusters with replication rate. Those set of clusters are high confidence set of poly(A) clusters. Poly(A) cluster .bed file have following format:


| Chromosome   |   Start |     End |   polyA_site |   count | Strand   | Feature         | gene_id            |      tpm |   gene_count |    usage |   fracA | signal         |   annotated_site |
|:-------------|--------:|--------:|-------------:|--------:|:---------|:----------------|:-------------------|---------:|-------------:|---------:|--------:|:---------------|-----------------:|
| chr17        | 2681887 | 2681912 |      2681907 |     249 | +        | three_prime_utr | ENSG00000007168.14 |  7978.72 |          780 | 0.319231 |       4 | 2681885@AATAAA |          2685608 |
| chr17        | 2684498 | 2684510 |      2684505 |      58 | +        | three_prime_utr | ENSG00000007168.14 |  1858.5  |          780 | 0.074359 |       4 | 2684477@AATAAA |          2685608 |
| chr17        | 2685607 | 2685616 |      2685614 |     473 | +        | three_prime_utr | ENSG00000007168.14 | 15156.4  |          780 | 0.60641  |       3 | 2685562@ATTAAA |          2685608 |
| chr17        | 3661532 | 3661541 |      3661536 |     110 | +        | three_prime_utr | ENSG00000040531.16 |  3524.74 |          110 | 1        |       1 | 3661514@ATTAAA |          3663103 |
| chr17        | 2059842 | 2059845 |      2059843 |     145 | -        | three_prime_utr | ENSG00000070366.14 |  4646.24 |          145 | 1        |       2 | 2059861@AATAAA |          2059842 |
...

This cluster .bed file and all other cluster .bed files have following columns:


- **Chromosome:** Chromosome of the poly(A) cluster.
- **Start:** Start position of poly(A) cluster.
- **End:** End position of poly(A) cluster.
- **polyA_site:** Exact poly(A) site (peak) of the cluster.
- **count:** number of reads supporting to cluster (ending in cluster).
- **Strand:** Strand of the poly(A) cluster.
- **Feature:** Genomics feature overlapping with the cluster (obtained from GTF file).
- **gene_id:** The gene containing the poly(A) clusters.
- **tpm:** TPM of the cluster calculated by $count / sum(count) * 1,000,000$.
- **gene_count:** Total number reads in all the clusters of this gene calculated by $sum(count\_i)$ where $i \in gene$.
- **usage:** Percentage use of specific poly(A) clusters of the gene calculated by $count / gene\_count$
- **fracA:** Number of A bp following the poly(A) cluster peak indicationg possible internal priming.
- **Signal:** Poly(A) signal of the cluster have format of start position of the cluster and signal sequence separated by @.
- **annotated_site:** End position of the 3' UTR based on the GTF if poly(A) cluster located in 3' UTR.

Poly(A) clusters can be read by following code as dataframe:

```
from lapa import read_polyA_cluster

df = read_polyA_cluster('your_output_dir/polyA_clusters.bed')
```

`raw_polyA_clusters.bed` contains the all the poly(A) clusters detected by LAPA but not filtered for replication.


`counts` is directory containing read end bigwig files. Each bigwig file contains number of reads ends per position indicating possosible poly(A) sites. This directory contains one bigwig file for each strand and not filtered so representing row data. There is pair of bigwig files per sample. The file starting all prefix contains counts from all the samples where counts are aggreated into one bigwig file.


`coverage` is directory containing bigwig files for coverage. Each file indicates for coverage of non-zero values in counts file. So the file format is sparse, contains values only for positions where at least 1 read is ending, and remaining positions are zero despite coverage can be non-zero. Sparse file format used to limit file and computational efficiency. The file starting all prefix contains counts from all the samples where counts are aggreated into one bigwig file.


`coverage` is directory containing ratio of counts to coverage ($ count / coverage $). This ratio indicates percentage reads are ending at a position given coverage. If the ratio close to one, the site is definitive poly(A) site give there is high coverage. If ratio is close to 0 then reads ending could be at the position by chance. Based on the default parameters LAPA (`cluster_ratio_cutoff`) only initialize cluster if ratio > 5% at a position. The file starting all prefix contains counts from all the samples where counts are aggreated into one bigwig file.


`dataset` is directory containing poly(A) cluster .bed files per dataset. Those files are filtered for replication rate using samples in for the dataset, then replicated clusters from all the samples aggregate into .bed file for the cluster.


`raw_sample` is directory containing poly(A) cluster .bed files per sample where files are not filtered for replication. 

`sample` is directory containing poly(A) cluster .bed files per sample where files are filtered for replication. 


`logs` is directory containing logs of LAPA. `final_stats.log` contains statistics about poly(A) clusters after program finished. `progress.log` provide inside about the progress of program run. `warnings.log` file contains possible warning encounter during the run time if there is any.


## TSS clusters

Output directory of TSS LAPA (`lapa_tss --output_dir your_output_dir`) looks like:

```
your_output_dir/
├── tss_clusters.bed
├── raw_tss_clusters.bed
├── counts
│   ├── all_polyA_counts_neg.bw
│   ├── all_polyA_counts_pos.bw
│   ├── {sample}_polyA_counts_neg.bw
│   ├── {sample}_polyA_counts_pos.bw
├── coverage
│   ├── all_polyA_coverage_neg.bw
│   ├── all_polyA_coverage_pos.bw
│   ├── {sample}_polyA_coverage_neg.bw
│   ├── {sample}_polyA_coverage_pos.bw
├── ratio
│   ├── all_polyA_ratio_neg.bw
│   ├── all_polyA_ratio_pos.bw
│   ├── {sample}_polyA_ratio_neg.bw
│   ├── {sample}_polyA_ratio_pos.bw
├── dataset
│   └── {dataset}.bed
├── raw_sample
│   └── {sample}.bed
├── sample
│   └── {sample}.bed
├── logs
│   ├── final_stats.log
└── ├── progress.log
    └── warnings.log
```

where `tss_clusters.bed` is the main output of LAPA and contains TSS clusters with replication rate. Those set of clusters are high confidence set of TSS clusters. TSS cluster .bed file have following format:


| Chromosome   |    Start |      End |   tss_site |   count | Strand   | Feature        | gene_id            |     tpm |   gene_count |     usage |   annotated_site |
|:-------------|---------:|---------:|-----------:|--------:|:---------|:---------------|:-------------------|--------:|-------------:|----------:|-----------------:|
| chr17        | 38870035 | 38870063 |   38870060 |     201 | +        | five_prime_utr | ENSG00000002834.18 | 2220.6  |          252 | 0.797619  |         38869858 |
| chr17        | 38890455 | 38890456 |   38890456 |      24 | +        | exon           | ENSG00000002834.18 |  265.15 |          252 | 0.0952381 |               -1 |
| chr17        | 38918997 | 38918998 |   38918998 |      27 | +        | exon           | ENSG00000002834.18 |  298.29 |          252 | 0.107143  |               -1 |
| chr17        | 48107566 | 48107599 |   48107576 |      13 | +        | five_prime_utr | ENSG00000002919.15 |  143.62 |           23 | 0.565217  |         48107548 |
| chr17        | 48107764 | 48107785 |   48107785 |      10 | +        | five_prime_utr | ENSG00000002919.15 |  110.48 |           23 | 0.434783  |         48107548 |
...

This cluster .bed file and all other cluster .bed files have following columns:


- **Chromosome:** Chromosome of the poly(A) cluster.
- **Start:** Start position of poly(A) cluster.
- **End:** End position of poly(A) cluster.
- **tss_site:** Exact poly(A) site (peak) of the cluster.
- **count:** number of reads supporting to cluster (ending in cluster).
- **Strand:** Strand of the poly(A) cluster.
- **Feature:** Genomics feature overlapping with the cluster (obtained from GTF file).
- **gene_id:** The gene containing the poly(A) clusters.
- **tpm:** TPM of the cluster calculated by $count / sum(count) * 1,000,000$.
- **gene_count:** Total number reads in all the clusters of this gene calculated by $sum(count\_i)$ where $i \in gene$.
- **usage:** Percentage use of specific poly(A) clusters of the gene calculated by $count / gene\_count$
- **annotated_site:** End position of the 5' UTR based on the GTF if poly(A) cluster located in 5' UTR.

`raw_tss_clusters.bed` contains the all the TSS clusters detected by LAPA but not filtered for replication.

For the details of the other files see the documentation of Poly(A) clusters. File structure and contain of the files are same with TSS clusters.
