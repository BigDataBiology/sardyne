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
    '''Extract gene positions from a Prodigal FASTA file'''
    data = []
    for line in open(f):
        if match := re.match(PRODIGAL_FASTA_HEADER_PAT, line):
            gene_id, start, end, is_reverse = match.groups()
            start = int(start)
            end = int(end)
            is_reverse = is_reverse == '-'
            length = abs(end - start)
            data.append((gene_id, start, end, length, is_reverse))
    return pl.DataFrame(data, schema={
                    'gene_id': pl.String,
                    'start': pl.Int32,
                    'end': pl.Int32,
                    'length': pl.Int32,
                    'is_reverse': pl.Boolean})


def load_diamond_outputs(diamond_outs):
    '''Load DIAMOND output files

    Parameters
    ----------
    diamond_outs: list of paths to DIAMOND output files

    Returns
    -------
    polars.DataFrame
    '''
    diamond_out = pl.concat([
        pl.read_csv(f, separator='\t', has_header=False)
            for f in diamond_outs])
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
    return diamond_out


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
    view_muts = [0, 100, 1000]
    diamond_out = load_diamond_outputs(glob(f'outputs/checkm2_{tag}_simulation/diamond_output/DIAMOND_RESULTS*.tsv'))
    diamond_out = diamond_out.with_columns([
        diamond_out['genome'].map_elements(
                                lambda g: int(g.split('_')[1].split('.')[0]),
                                return_dtype=pl.Int32
                                ).alias('nr_mutations')
        ])
    data = []
    for mut in view_muts:
        diamond_out_mut = diamond_out.filter(pl.col('nr_mutations') == mut)
        f = f'outputs/checkm2_{tag}_simulation/protein_files/mutated_{mut:06}.fna.faa'
        orfs = get_gene_positions(f)
        diamond_out_mut = diamond_out_mut.join(orfs, right_on='gene_id', left_on='orf')
        data.append(diamond_out_mut[['nr_mutations', 'KO', 'length']])
    data = pl.concat(data)
    data = data.filter(pl.col('KO').is_in(aa_sizes_by_ko.keys()))

    zscores = []
    for ix in range(len(data)):
        nr_muts, ko, s = data.row(ix)
        db_size = np.array(aa_sizes_by_ko[ko])*3

        zscores.append((nr_muts,
                        ko,
                        (s - db_size.mean())/db_size.std(),
                        ))

    zscores = pl.DataFrame(zscores, schema={
        'nr_muts': pl.Int32,
        'KO': pl.String,
        'z': pl.Float32,
        })
    fig, ax = plt.subplots()
    sns.boxplot(data=zscores,
                x='nr_muts',
                y='z',
                boxprops=dict(alpha=.4, facecolor='w'),
                width=.2,
                ax=ax,
                )
    ax.set_title(f'z-scores for KO lengths ({tag})')
    fig.savefig(f'plots/length_by_ko_checkm2_zscore_{tag}.png')

    fig,ax = plt.subplots()
    X = np.linspace(-12,12,1000)
    for lim,lab in [(0, 'WT'), (100, '100'), (1000, '1k')]:
        sel = zscores.filter(pl.col('nr_muts') == lim)
        ax.plot(X, [(sel['z'] < x).mean() for x in X], label=lab)
    ax.legend()
    ax.set_xlabel('KO length z-score')
    ax.set_ylabel('Fraction of genes (cumm)')
    fig.tight_layout()
    fig.savefig(f'plots/ko_zscores_{tag}.png')

    zscores_thresh = pl.concat([
                zscores.with_columns([
                        (zscores['z'] < lim).alias('bel_lim')
                        ]) \
                    .group_by('nr_muts') \
                    .sum()[['nr_muts', 'bel_lim']] \
                    .with_columns([pl.lit(lim).alias('lim')])
                for lim in [-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12]])
    zscores_thresh = zscores_thresh.pivot(columns=['nr_muts'], values=['bel_lim'], index=['lim'])
    print(tag)
    print(zscores_thresh)
