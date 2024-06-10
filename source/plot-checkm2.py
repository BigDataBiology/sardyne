from os import makedirs
import numpy as np
import polars as pl
from matplotlib import pyplot as plt
from matplotlib import cm

from glob import glob
import jug

colors = cm.Dark2.colors

_, jugspace = jug.init('simulate.py', 'simulate.jugdata')
gene_sizes = jug.value(jugspace['gene_sizes'])
nr_muts = jugspace['nr_muts']

makedirs('plots', exist_ok=True)
for tag in gene_sizes:
    q = pl.read_csv(
            f'./outputs/checkm2_{tag}_simulation/quality_report.tsv',
            separator='\t')
    q = q.with_columns([
        q['Name'] \
                .map_elements(lambda s:int(s.split('.')[0].split('_')[-1]), return_dtype=pl.Int32) \
                .rename("nr_mut")
        ])
    q = q.sort('nr_mut')

    avg_gene_size = np.array([
        gs.mean() for gs in gene_sizes[tag]])

    fig,ax1 = plt.subplots()
    ax1.plot(q['nr_mut'], q['Completeness'],'o-', ms=1, c=colors[0], label='Completeness (checkM2)')
    ax1.set_xlabel('Nr mutations')
    ax1.set_ylabel('Completeness')

    ax2 = ax1.twinx()
    ax2.plot(nr_muts, avg_gene_size, c=colors[1], label='Avg. gene size')
    ax2.set_ylabel('Avg. gene size')
    fig.legend()

    fig.savefig(f'plots/check2_{tag}_mutations_50k.png')

    ax1.set_xlim(0, 5_000)
    fig.savefig(f'plots/check2_{tag}_mutations_5k.png')

    data = pl.concat([ pl.read_csv(f, separator='\t', has_header=False)
        for f in
            glob(f'outputs/checkm2_{tag}_simulation/diamond_output/DIAMOND_RESULTS_*.tsv')])
    data.columns = ['query',
                    'subject',
                    'identity',
                    'alignment_length',
                    'mismatches',
                    'gap_opens',
                    'q_start',
                    'q_end',
                    's_start',
                    's_end',
                    'evalue',
                    'bit_score']

    data = data.with_columns([
        data['query'].map_elements(lambda q: q.split('Ω')[0], return_dtype=str).alias('genome'),
        data['subject'].map_elements(lambda s: s.split('~')[1], return_dtype=str).alias('KO')
        ])

    genome2kos = data.group_by('genome').agg(pl.col('KO').map_elements(set, return_dtype=pl.Object))
    S = genome2kos.to_dict()
    genome2kos = dict(zip(S['genome'], S['KO']))
    nr_muts2kos = {int(k.split('_')[1].split('.')[0]): v for k, v in genome2kos.items()}
    assert nr_muts == sorted(nr_muts2kos.keys())

    unmutatate_kos = nr_muts2kos[0]
    false_positives = [len(nr_muts2kos[n] - unmutatate_kos) for n in nr_muts]
    false_negatives = [len(unmutatate_kos - nr_muts2kos[n]) for n in nr_muts]

    fig, ax = plt.subplots()
    ax.plot(nr_muts, false_negatives, label='False negatives', c=colors[0])
    ax.set_xlabel('Nr mutations')
    ax.set_ylabel('Nr missing KOs')

    ax2 = ax.twinx()
    ax2.plot(nr_muts, false_positives, label='False positives', c=colors[1])
    ax2.set_ylabel('Nr false positive KOs')
    fig.legend()

    fig.savefig(f'plots/checkm2_{tag}_false_positives_negatives_50k.png')
    ax.set_xlim(0, 5_000)
    fig.savefig(f'plots/checkm2_{tag}_false_positives_negatives_5k.png')

