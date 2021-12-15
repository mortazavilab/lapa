from matplotlib import pyplot as plt
from lapa.count import TailTesCounter

ttc = TailTesCounter(snakemake.input[0], min_tail_len=1)
ttc.plot_tail_len_dist()
plt.savefig(snakemake.output[0])

