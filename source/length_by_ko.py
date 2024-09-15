from os import makedirs
import gzip
import polars as pl
import numpy as np
from collections import defaultdict
import seaborn as sns
from matplotlib import pyplot as plt
from fasta import fasta_iter
import re
from glob import glob
import jug

PRODIGAL_FASTA_HEADER_PAT = r'^>(\S+) # (\d+) # (\d+) # (-?)1 # '
MIN_UNIGENES_PER_KO = 100

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

with gzip.open('./outputs/GMGC10.complete.annotated.tsv.gz', 'rb') as f:
    gmgc = pl.read_csv(f, separator='\t')

gmgc = gmgc.with_columns([
    gmgc['subject'].str.split('~').list[1].alias('KO'),
    ])
kos = set(gmgc.group_by('KO').len().filter(pl.col('len') > MIN_UNIGENES_PER_KO)['KO'].to_list())
pre_len = len(gmgc)
gmgc = gmgc.filter(pl.col('KO').is_in(kos))
post_len = len(gmgc)
print(f'Removed {pre_len - post_len:,} unigenes (out of {pre_len:,}; {(pre_len - post_len)/pre_len:.2%})')

aa_sizes_by_ko = defaultdict(list)
for ko, unigene_len in gmgc[['KO', 'gmgc_unigene_len_aa']].iter_rows():
    aa_sizes_by_ko[ko].append(unigene_len)

aa_sizes_by_ko = {k:np.array(v) for k,v in aa_sizes_by_ko.items()}
for v in aa_sizes_by_ko.values():
    v.sort()

nr_seqs_by_ko = dict(gmgc.group_by('KO').len().iter_rows())

for tag,_ in jugspace['INPUT_DATA']:
    f_wt = f'outputs/checkm2_{tag}_simulation/protein_files/mutated_000000.fna.faa'
    f_100 = f'outputs/checkm2_{tag}_simulation/protein_files/mutated_000100.fna.faa'
    f_1k = f'outputs/checkm2_{tag}_simulation/protein_files/mutated_001000.fna.faa'

    diamond_out = pl.concat([ pl.read_csv(f, separator='\t', has_header=False)
        for f in
            glob(f'outputs/checkm2_{tag}_simulation/diamond_output/DIAMOND_RESULTS*.tsv')])
    diamond_out.columns = ['query',
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

    diamond_out = diamond_out.with_columns([
        diamond_out['query'].str.split('Ω').list[0].alias('genome'),
        diamond_out['query'].str.split('Ω').list[1].alias('orf'),
        diamond_out['subject'].str.split('~').list[1].alias('KO')
        ])

    wt = diamond_out.filter(pl.col('genome') == 'mutated_000000.fna')
    mut100 = diamond_out.filter(pl.col('genome') == 'mutated_000100.fna')
    mut1k = diamond_out.filter(pl.col('genome') == 'mutated_001000.fna')

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

    zscores = []
    for ix in range(len(data)):
        ko, s_wt, s_100, s_1k = data.row(ix)
        if ko not in aa_sizes_by_ko:
            continue
        db_size = np.array(aa_sizes_by_ko[ko])*3

        zscores.append((ko,
            (s_wt - db_size.mean())/db_size.std(),
            (s_100 - db_size.mean())/db_size.std(),
            (s_1k - db_size.mean())/db_size.std(),
                        ))

    zscores = pl.DataFrame(zscores, schema={'KO': pl.String,
                                        'z_wt': pl.Float32,
                                        'z_100': pl.Float32,
                                        'z_1k': pl.Float32})
    fig, ax = plt.subplots()
    sns.boxplot(data=zscores.melt(id_vars=['KO']),
                x='variable',
                y='value',
                boxprops=dict(alpha=.4, facecolor='w'),
                width=.2,
                ax=ax,
                )
    ax.set_title(f'z-scores for KO lengths ({tag})')
    fig.savefig(f'plots/length_by_ko_checkm2_zscore_{tag}.png')

    fig,ax = plt.subplots()
    X = np.linspace(-12,12,1000)
    ax.plot(X, [(zscores['z_wt'] < x).mean() for x in X], label='WT')
    ax.plot(X, [(zscores['z_100'] < x).mean() for x in X], label='100')
    ax.plot(X, [(zscores['z_1k'] < x).mean() for x in X], label='1k')
    ax.legend()
    ax.set_xlabel('KO length z-score')
    ax.set_ylabel('Fraction of genes (cumm)')
    fig.tight_layout()
    fig.savefig(f'plots/ko_zscores_{tag}.png')

    zscores_thresh = pl.concat([
                (zscores.select(pl.col(pl.Float32)) < lim).sum().with_columns(lim=lim)
                    for lim in [-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12]])
    print(tag)
    print(zscores_thresh)
