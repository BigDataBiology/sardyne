from os import makedirs
import gzip
import polars as pl
import numpy as np
from collections import defaultdict
import seaborn as sns
from matplotlib import pyplot as plt
import re
from glob import glob
import jug

MIN_UNIGENES_PER_KO = 100

_, jugspace = jug.init('simulate.py', 'simulate.jugdata')
makedirs('plots', exist_ok=True)

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
    diamond_out = diamond_out.with_columns([
        diamond_out['genome'].map_elements(
                                lambda g: int(g.split('_')[1].split('.')[0]),
                                return_dtype=pl.Int32
                                ).alias('nr_mutations')
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

ko_sizes = pl.DataFrame(
        [ # Multiply by 3 to get nucleotide length
            [k, v.mean() * 3, (v*3).std()]
            for k,v in aa_sizes_by_ko.items()],
        orient='row',
        schema={
            'KO': pl.String,
            'mean': pl.Float32,
            'std': pl.Float32,
            }
        )

view_muts = [0, 100, 1000]
zscores_outs = []
esgs_fig,esgs_ax = plt.subplots()

for tag,_ in jugspace['INPUT_DATA']:
    diamond_out = load_diamond_outputs(glob(f'outputs/checkm2_{tag}_simulation/diamond_output/DIAMOND_RESULTS*.tsv'))
    all_muts = sorted(set(diamond_out['nr_mutations']))
    with gzip.open(f'outputs/checkm2_{tag}_simulation/all_gene_positions.tsv.gz', 'rb') as f:
        all_gene_positions = pl.read_csv(f, separator='\t')

    data = diamond_out.join(all_gene_positions, right_on='gene_id', left_on='orf')
    data = data.filter(pl.col('KO').is_in(aa_sizes_by_ko.keys()))

    zscores = data.join(ko_sizes, left_on='KO', right_on='KO').with_columns(
        z=(pl.col('length')- pl.col('mean'))/pl.col('std')
        )[['nr_mutations', 'KO', 'z']]

    esgs = []
    for mut in all_muts:
        sel = zscores.filter((pl.col('nr_mutations') == mut) & (pl.col('z') < -4)).select(pl.col('z').sum())
        esgs.append(sel.item())

    [ell] = esgs_ax.plot(all_muts, esgs, label=tag)
    esgs_ax.plot(all_muts, [esgs[0] for _ in all_muts], label=None, linestyle='--', color=ell.get_color())
    esgs_ax.set_xlabel('Number of mutations')
    esgs_ax.set_ylabel('Sum of z-scores < -4')

    view_muts = [0, 100, 1000]
    zscores = zscores.filter(pl.col('nr_mutations').is_in(view_muts))

    fig, ax = plt.subplots()
    sns.boxplot(data=zscores,
                x='nr_mutations',
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
        sel = zscores.filter(pl.col('nr_mutations') == lim)
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
                    .group_by('nr_mutations') \
                    .sum()[['nr_mutations', 'bel_lim']] \
                    .with_columns([pl.lit(lim).alias('lim')])
                for lim in [-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12]])
    zscores_thresh = zscores_thresh.pivot(on=['nr_mutations'], values=['bel_lim'], index=['lim'])
    zscores_outs.append(zscores_thresh.with_columns(tag=pl.lit(tag)))

esgs_ax.legend()
sns.despine(esgs_fig, trim=True)
esgs_fig.tight_layout()
esgs_fig.savefig('plots/esgs_m4.svg')
esgs_fig.savefig('plots/esgs_m4.png')
