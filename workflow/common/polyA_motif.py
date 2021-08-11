import matplotlib.pyplot as plt
from Bio import motifs
from Bio.Seq import Seq
from lapa.utils.common import polyA_signal_seqs


instances = [Seq(i) for i in polyA_signal_seqs]
motif = motifs.create(instances)

print(motif.degenerate_consensus)
plt.figure(figsize=(10, 5), dpi=250)
motif.weblogo(snakemake.output['motif_logo'], format="pdf")
