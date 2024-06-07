from matplotlib import pyplot as plt
import numpy as np
import polars as pl
from matplotlib import cm
colors = cm.Dark2.colors
import jug

_, jugspace = jug.init('simulate.py', 'simulate.jugdata')
gene_sizes = jug.value(jugspace['gene_sizes'])
nr_muts = jugspace['nr_muts']

q = pl.read_csv(
        './outputs/checkm2_ecoli_k12_simulation/quality_report.tsv',
        separator='\t')
q = q.with_columns([
    q['Name'] \
            .map_elements(lambda s:int(s.split('.')[0].split('_')[-1]), return_dtype=pl.Int32) \
            .rename("nr_mut")
    ])
q = q.sort('nr_mut')

avg_gene_size = np.array([
    gs.mean() for gs in gene_sizes['ecoli_k12']])

fig,ax1 = plt.subplots()
ax1.plot(q['nr_mut'], q['Completeness'],'o-', ms=1, c=colors[0], label='Completeness (checkM2)')
ax1.set_xlabel('Nr mutations')
ax1.set_ylabel('Completeness')

ax2 = ax1.twinx()
ax2.plot(nr_muts, avg_gene_size, c=colors[1], label='Avg. gene size')
ax2.set_ylabel('Avg. gene size')
fig.legend()

fig.savefig('check2_mutations_50k.png')

ax1.set_xlim(0, 5_000)
fig.savefig('check2_mutations_5k.png')
