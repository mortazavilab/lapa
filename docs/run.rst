.. _run-lapa:
Running LAPA with CLI and python
===================================


.. _run-lapa-polya:
Poly(A) cluster calling
------------------------

.. code-block:: bash

   lapa --alignment {bam} \
	--fasta {fasta} \
	--annotation {gtf} \
	--chrom_sizes {chrom_sizes} \
	--output_dir {output}

Output of LAPA (``--output_dir``) is directory and content is explained in :ref:`output format section of documentation<output-format>`.

Fasta is standart genome sequence file can be downloaded from Ensembl or GENCODE. Based on the FASTA file `chrom_sizes` file can be generated with the following command with faidx:


.. code-block:: bash

   pip install pyfaidx
   faidx {fasta} -i chromsizes > {chrom_sizes}

GTF similarly can be downloaded from (Ensembl or GENCODE). Ensembl separately annotates UTR features of GTF file as **three_prime_utr** and **five_prime_utr** so directly used as input to LAPA. Yet GENCODE does not distinguish 5' UTR or 3' UTR and annotated all as just **UTR**; thus, UTR entry of GENCODE needs to be corrected before passing it to LAPA with the following command with a tool called `gencode_utr_fix`:

.. code-block:: bash

   pip install cython
   pip install -e git+https://github.com/MuhammedHasan/gencode_utr_fix.git#egg=gencode_utr_fix
   gencode_utr_fix --input_gtf {gencode_gtf} --output_gtf {gencode_utr_fixed_gtf}

LAPA can process single or multiple bam files as input. For example, single bam file:

.. code-block:: bash

   lapa --alignment single.bam \
	--fasta {fasta} \
	--annotation {gtf} \
	--chrom_sizes {chrom_sizes} \
	--output_dir {output}

.. warning::
   Running LAPA jointly with multiple samples (bam files) is highly recommended over running it on individual samples (bam files) separately.

Multiple bam files can be given as input with comma separation:

.. code-block:: bash

   lapa --alignment sample1.bam,sample2.bam,sample3.bam \
	--fasta {fasta} \
	--annotation {gtf} \
	--chrom_sizes {chrom_sizes} \
	--output_dir {output}

If multiple bam files are given as input, LAPA filters poly(A) clusters based on the replication rate between samples (:ref:`See the documentation for CLI<cli-arguments-lapa>`). If comma separated bam files are provided with the command above, LAPA assumes all the files are experimental replicates of each other.

For multiple samples, CSV file can also be provided as alignment input then replication rate will be calculated based on the dataset column:

.. list-table:: sample.csv
   :widths: 10 10 10
   :header-rows: 1

   * - sample
     - dataset
     - path
   * - biosample1_rep1
     - tissue1
     - your_bam_dir/bam_sample1_rep1_path.bam
   * - biosample1_rep2
     - tissue1
     - your_bam_dir/bam_sample1_rep2_path.bam
   * - biosample2_rep1
     - tissue2
     - your_bam_dir/bam_sample2_rep1_path.bam
   * - biosample2_rep2
     - tissue2
     - your_bam_dir/bam_sample2_rep1_path.bam

Rows of the file is separated by comma:

.. code-block::

   sample,dataset,path
   biosample1_rep1,tissue1,bam_sample1_rep1_path.bam
   biosample1_rep2,tissue1,bam_sample1_rep2_path.bam
   biosample2_rep1,tissue2,bam_sample2_rep1_path.bam
   biosample2_rep2,tissue2,bam_sample2_rep2_path.bam

Then LAPA takes sample.csv as alignment input:

.. code-block:: bash

   lapa --alignment samples.csv \
	--fasta {fasta} \
	--annotation {gtf} \
	--chrom_sizes {chrom_sizes} \
	--output_dir {output}

.. note::
   The replication rate is the most important parameter of this package because it determines the final set of poly(A) clusters by filtering samples below a specific replication rate. If you have experimental replicates, pass them in the sample file as a **dataset**. By default, LAPA assumes datasets are experimental replicates of each other and applies a highly conservative 0.95 replication rate cutoff. If experimental replicates are biological replicates with batch effect .etc, the default replication rate could be stringent. For a higher recall setting, the replication rate could be set to a lower cutoff such as (`--min_replication_rate 0.75`). Moreover, the default behavior of replication rate accepts clusters as replicated if there is at least 1 read 2 out of N samples. If you have many ~10-100 replicates, you may set `replication_num_sample` to a higher number just of 5 or 10% of samples. See the quality control section to visualize and choose the better replication rate. (:ref:`See the documentation for CLI<cli-arguments-lapa>`)


End counting and tail counting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



   
Calling python API:
~~~~~~~~~~~~~~~~~~~

LAPA also provides high level python api to call poly(A) clusters:

.. code-block:: python

    from lapa import lapa

    lapa(alignment='sample1.bam,sample2.bam,sample3.bam',
	 fasta='hg38.fasta',
	 annotation='hg38.gtf',
	 chrom_sizes='hg38.chrom_sizes',
	 output_dir='output_dir')

See the API reference for other options.


.. _run_lapa_tss:
TSS cluster calling
--------------------

Options for TSS cluster calling for similar so see :ref:`See the documentation for CLI<run-lapa-polya>`)


.. code-block:: bash

   lapa_tss --alignment sample1.bam,sample2.bam,sample3.bam \
	    --fasta {fasta} \
	    --annotation {gtf} \
	    --chrom_sizes {chrom_sizes} \
	    --output_dir {output}

Similarly, you can run ``lapa_tss`` on multiple samples
	
.. code-block:: bash

   lapa_tss --alignment samples.csv \
	    --fasta {fasta} \
	    --annotation {gtf} \
	    --chrom_sizes {chrom_sizes} \
	    --output_dir {output}

See (:ref:`See the documentation for CLI<cli-arguments-lapa-tss>`) further options of LAPA.

Calling python API:
~~~~~~~~~~~~~~~~~~~

LAPA also provides high level python api to call TSS clusters:

.. code-block:: python

    from lapa import lapa_tss

    lapa_tss(alignment='sample1.bam,sample2.bam,sample3.bam',
             fasta='hg38.fasta',
	     annotation='hg38.gtf',
	     chrom_sizes='hg38.chrom_sizes',
	     output_dir='output_dir')

See the API reference for other options.


Linking reads
---------------

After performing TSS and poly(A) cluster calling, linking reads can be calculated from long-reads. Linking reads are the reads start in TSS cluster and ends in poly(A) cluster. 

.. code-block:: bash

   lapa_link_tss_to_tes --alignment sample_rep1.bam \
		        --lapa_dir {lapa_dir} \
			--lapa_tss_dir {lapa_tss_dir} \
			--output {linking_reads}.csv


Dataset specific linking can be obtained with:

.. code-block:: bash

   lapa_link_tss_to_tes --alignment sample_rep1.bam \
		        --lapa_dir {lapa_dir} \
			--lapa_tss_dir {lapa_tss_dir} \
			--output {linking_reads}.csv \
			--dataset {dataset_name}

where {dataset_name} defined in lapa_dir and lapa_tss_dir.

See (:ref:`See the documentation for CLI<cli-arguments-lapa-link_tss-to-tes>`) further details of arguments.

Calling python API:
~~~~~~~~~~~~~~~~~~~

LAPA also provides high level python api to call TSS clusters:

.. code-block:: python

    from lapa.link import link_tss_to_tes

    df_links = link_tss_to_tes(alignment='sample1.bam',
		               lapa_dir='hg38.fasta',
			       lapa_tss_dir='hg38.gtf',
			       dataset='{dataset_name}')

See the API reference for other options.


Correction of TALON gtf and abundance
---------------------------------------

LAPA creates GTF file with tss/poly(A) cluster support based on the linking reads and using splice chain of TALON:

.. code-block:: bash

   lapa_correct_talon --links {links_reads}.csv \
		      --read_annot {talon_read_annot} \
		      --gtf_input {talon_gtf} \
		      --abundance_input {talon_abundance} \
		      --gtf_output {corrected_gtf} \
		      --abundance_output {corrected_abundance}

where ``abundance``, ``gtf``, ``read_annot`` files are defined with TALON.

See (:ref:`See the documentation for CLI<cli-arguments-lapa-correct-talon>`) further details of arguments.

LAPA easly can be adapted to start and end correction of the tools beyond TALON. See the linking read api for details and issue/pull request if you want to integrate LAPA to your long-read tools.

Calling python API:
~~~~~~~~~~~~~~~~~~~

LAPA also provides high level python api to call TSS clusters:

.. code-block:: python

    from lapa.correction import correct_talon

    df_links = link_tss_to_tes(link_path='linking_reads.csv'
		               read_annot_path='talon_read_annot.tsv',
			       gtf_input='talon_input.gtf',
			       gtf_output='talon_corected.gtf',
			       abundance_path='talon_abundance.csv',
			       abundance_output='talon_abundance_corrected.csv')
