from os import makedirs
import polars as pl
import numpy as np
from collections import defaultdict
import seaborn as sns
from matplotlib import pyplot as plt
from fasta import fasta_iter
import re
from glob import glob
import jug

from fasta import fasta_iter

PRODIGAL_FASTA_HEADER_PAT = r'^>(\S+) # (\d+) # (\d+) # (-?)1 # '

_, jugspace = jug.init('simulate.py', 'simulate.jugdata')
makedirs('plots', exist_ok=True)

def get_gene_positions(f):
    data = []
    for line in open(f):
        if match := re.match(PRODIGAL_FASTA_HEADER_PAT, line):
            gene_id, start, end, is_reverse = match.groups()
            start = int(start)
            end = int(end)
            is_reverse = is_reverse == '-'
            data.append((gene_id, start, end, is_reverse))
    return pl.DataFrame(data, schema={
                    'gene_id': pl.String,
                    'start': pl.Int32,
                    'end': pl.Int32,
                    'is_reverse': pl.Boolean})

aa_sizes_by_ko = defaultdict(list)
for h,seq in fasta_iter('../data/data/uniref100.KO.faa'):
    _, ko = h.split('~')
    aa_sizes_by_ko[ko].append(len(seq))
nr_seqs_by_ko = {k:len(orfs) for k,orfs in aa_sizes_by_ko.items()}


for tag,_ in jugspace['INPUT_DATA']:
    f_wt = f'outputs/checkm2_{tag}_simulation/protein_files/mutated_000000.fna.faa'
    f_100 = f'outputs/checkm2_{tag}_simulation/protein_files/mutated_000100.fna.faa'
    f_1k = f'outputs/checkm2_{tag}_simulation/protein_files/mutated_001000.fna.faa'

    data = pl.concat([ pl.read_csv(f, separator='\t', has_header=False)
        for f in
            glob(f'outputs/checkm2_{tag}_simulation/diamond_output/DIAMOND_RESULTS*.tsv')])
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
        data['query'].map_elements(lambda q: q.split('Ω')[1], return_dtype=str).alias('orf'),
        data['subject'].map_elements(lambda s: s.split('~')[1], return_dtype=str).alias('KO')
        ])

    wt = data.filter(pl.col('genome') == 'mutated_000000.fna')
    mut100 = data.filter(pl.col('genome') == 'mutated_000100.fna')
    mut1k = data.filter(pl.col('genome') == 'mutated_001000.fna')

    wt_kos = set(wt['KO'].to_list())
    mut100_kos = set(mut100['KO'].to_list())
    mut1k_kos = set(mut1k['KO'].to_list())
    orfs_wt = get_gene_positions(f_wt)
    orfs_100 = get_gene_positions(f_100)
    orfs_1k = get_gene_positions(f_1k)
    data = []
    for k in (wt_kos & mut1k_kos & mut100_kos):
        wt_sel = wt.filter(pl.col('KO') == k)
        mut100_sel = mut100.filter(pl.col('KO') == k)
        mut1k_sel = mut1k.filter(pl.col('KO') == k)
        if len(wt_sel) == 1 and len(mut100_sel) == 1 and len(mut1k_sel) == 1:
            [wt_orf] = wt_sel['orf']
            [mut100_orf] = mut100_sel['orf']
            [mut1k_orf] = mut1k_sel['orf']
            _, start, end, _ = orfs_wt.filter(pl.col('gene_id') == wt_orf).row(0)
            _, start100, end100, _ = orfs_100.filter(pl.col('gene_id') == mut100_orf).row(0)
            _, start1k, end1k, _ = orfs_1k.filter(pl.col('gene_id') == mut1k_orf).row(0)
            data.append((k, end - start, end100 - start100, end1k - start1k))


    data = pl.DataFrame(data, schema={'KO': pl.String,
                                      'length_wt': pl.Int32,
                                      'length_100': pl.Int32,
                                      'length_1k': pl.Int32})
    diff = data.filter(pl.col('length_wt') != pl.col('length_1k')).sort('KO')

    fig, ax = plt.subplots()
    sns.boxplot(data=diff.melt(id_vars=['KO']),
                    x='variable',
                    y='value',
                    boxprops=dict(alpha=.4, facecolor='w'),
                    width=.2,
                    ax=ax,
                )

    for _,s_wt,s_100,s_1k in diff.iter_rows():
        color = 'k'
        if s_wt < s_1k:
            color = 'r'
        elif abs(s_wt - s_1k) > 100:
            color = 'g'
        elif abs(s_wt - s_1k) < 10:
            continue
        ax.plot([0,1,2], [s_wt, s_100, s_1k], '-o', c=color, lw=0.7, alpha=0.3)

    fig.savefig(f'plots/length_by_ko_checkm2_{tag}.png')

    data = data.with_columns([
        data['KO'].map_elements(nr_seqs_by_ko.get, return_dtype=pl.Int32).rename('nr_seqs')
        ])
    percs = []
    for ix in range(len(data)):
        ko, s_wt, s_100, s_1k, nr_seqs = data.row(ix)
        if nr_seqs < 100:
            continue
        db_size = np.array(aa_sizes_by_ko[ko])*3
        db_size.sort()
        percs.append((ko,
            np.searchsorted(db_size, s_wt)/len(db_size),
            np.searchsorted(db_size, s_100)/len(db_size),
            np.searchsorted(db_size, s_1k)/len(db_size),
                      ))

    percs = pl.DataFrame(percs, schema={'KO': pl.String,
                                        'perc_wt': pl.Float32,
                                        'perc_100': pl.Float32,
                                        'perc_1k': pl.Float32})
    fig, ax = plt.subplots()
    sns.boxplot(data=percs.melt(id_vars=['KO']),
                x='variable',
                y='value',
                boxprops=dict(alpha=.4, facecolor='w'),
                width=.2,
                ax=ax,
                )
    fig.savefig(f'plots/length_by_ko_checkm2_perc_{tag}.png')

    fig,ax = plt.subplots()
    X = np.linspace(0,1,1000)
    ax.plot(X, [(percs['perc_wt'] < x).mean() for x in X], label='WT')
    ax.plot(X, [(percs['perc_100'] < x).mean() for x in X], label='100')
    ax.plot(X, [(percs['perc_1k'] < x).mean() for x in X], label='1k')
    ax.plot(X, X, ':k', label='y=x')
    ax.legend()
    ax.set_xlabel('KO length percentile')
    ax.set_ylabel('Fraction of genes (cumm)')
    fig.tight_layout()
    fig.savefig(f'plots/ko_percentile_{tag}.png')

