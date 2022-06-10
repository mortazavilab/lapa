.. _output-format:

Output format
==============


Poly(A) clusters
-----------------

Output directory of LAPA (``lapa --output_dir your_output_dir``) looks like:

.. code-block::

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

**polyA_clusters.bed:** is the main output of LAPA and contains poly(A) clusters with replication rate. Those set of clusters are high confidence set of poly(A) clusters. Poly(A) cluster .bed file have following format:

.. mdinclude:: table_polya_cluster.md

Poly(A) clusters can be read by following code as dataframe:

.. code-block::
   
   from lapa import read_polyA_cluster
   
   df = read_polyA_cluster('your_output_dir/polyA_clusters.bed')


**raw_polyA_clusters.bed:** contains the all the poly(A) clusters detected by LAPA but not filtered for replication.


**counts:** is directory containing read end bigwig files. Each bigwig file contains number of reads ends per position indicating possosible poly(A) sites. This directory contains one bigwig file for each strand and not filtered so representing row data. There is pair of bigwig files per sample. The file starting all prefix contains counts from all the samples where counts are aggreated into one bigwig file.


**coverage:** is directory containing bigwig files for coverage. Each file indicates for coverage of non-zero values in counts file. So the file format is sparse, contains values only for positions where at least 1 read is ending, and remaining positions are zero despite coverage can be non-zero. Sparse file format used to limit file and computational efficiency. The file starting all prefix contains counts from all the samples where counts are aggreated into one bigwig file.


**coverage:** is directory containing ratio of counts to coverage ($ count / coverage $). This ratio indicates percentage reads are ending at a position given coverage. If the ratio close to one, the site is definitive poly(A) site give there is high coverage. If ratio is close to 0 then reads ending could be at the position by chance. Based on the default parameters LAPA (**cluster_ratio_cutoff**) only initialize cluster if ratio > 5% at a position. The file starting all prefix contains counts from all the samples where counts are aggreated into one bigwig file.


**dataset:** is directory containing poly(A) cluster .bed files per dataset. Those files are filtered for replication rate using samples in for the dataset, then replicated clusters from all the samples aggregate into .bed file for the cluster.


**raw_sample:** is directory containing poly(A) cluster .bed files per sample where files are not filtered for replication. 

**sample:** is directory containing poly(A) cluster .bed files per sample where files are filtered for replication. 


**logs:** is directory containing logs of LAPA. **final_stats.log** contains statistics about poly(A) clusters after program finished. **progress.log** provide inside about the progress of program run. **warnings.log** file contains possible warning encounter during the run time if there is any.


TSS clusters
-------------

Output directory of TSS LAPA (`lapa_tss --output_dir your_output_dir`) looks like:

.. code-block::

    your_output_dir/
    ├── tss_clusters.bed
    ├── raw_tss_clusters.bed
    ├── counts
    │   ├── all_tss_counts_neg.bw
    │   ├── all_tss_counts_pos.bw
    │   ├── {sample}_tss_counts_neg.bw
    │   ├── {sample}_tss_counts_pos.bw
    ├── coverage
    │   ├── all_polyA_coverage_neg.bw
    │   ├── all_polyA_coverage_pos.bw
    │   ├── {sample}_tss_coverage_neg.bw
    │   ├── {sample}_tss_coverage_pos.bw
    ├── ratio
    │   ├── all_polyA_ratio_neg.bw
    │   ├── all_polyA_ratio_pos.bw
    │   ├── {sample}_tss_ratio_neg.bw
    │   ├── {sample}_tss_ratio_pos.bw
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


where ``tss_clusters.bed`` is the main output of LAPA and contains TSS clusters with replication rate. Those set of clusters are high confidence set of TSS clusters. TSS cluster .bed file have following format:

.. mdinclude:: table_tss_cluster.md

Tss clusters can be read by following code as dataframe:
	       
.. code-block::
   
   from lapa import read_tss_cluster
   
   df = read_polyA_cluster('your_output_dir/tss_clusters.bed')

This cluster .bed file and all other cluster .bed files have following columns:


- **Chromosome:** Chromosome of the tss cluster.
- **Start:** Start position of tss cluster.
- **End:** End position of tss cluster.
- **tss_site:** Exact tss site (peak) of the cluster.
- **count:** number of reads supporting to cluster (ending in cluster).
- **Strand:** Strand of the tss cluster.
- **Feature:** Genomics feature overlapping with the cluster (obtained from GTF file).
- **gene_id:** The gene containing the tss clusters.
- **tpm:** TPM of the cluster calculated by $count / sum(count) * 1,000,000$.
- **gene_count:** Total number reads in all the clusters of this gene calculated by $sum(count\_i)$ where $i \in gene$.
- **usage:** Percentage use of specific tss clusters of the gene calculated by $count / gene\_count$
- **annotated_site:** End position of the 5' UTR based on the GTF if tss cluster located in 5' UTR.

`raw_tss_clusters.bed` contains the all the TSS clusters detected by LAPA but not filtered for replication.

For the details of the other files see [the documentation of Poly(A) clusters](#poly(A)-clusters). File structure and content of the files are same with TSS clusters.
