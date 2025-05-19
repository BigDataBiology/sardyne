from os import makedirs
import gzip
import polars as pl
import numpy as np
from collections import defaultdict
import seaborn as sns
from matplotlib import pyplot as plt
import jug
from eggnog import extract_og

PLOT_KO_LENGTHS_BOXPLOT = True
PLOT_KO_LENGTHS_CUMMDIST = True
MIN_UNIGENES_PER_OG = 100

_, jugspace = jug.init('simulate.py', 'simulate.jugdata')
makedirs('plots/emapper_ref/', exist_ok=True)

with gzip.open('./outputs/GMGC10.emapper2.annotations.complete.tsv.gz', 'rb') as f:
    emapper = pl.read_csv(f, separator='\t')

ogs = set(
        emapper
            .group_by('eggNOG_OG1')
            .len()
            .filter(pl.col('len') > MIN_UNIGENES_PER_OG)['eggNOG_OG1']
            .to_list()
    )
pre_len = len(emapper)
emapper = emapper.filter(pl.col('eggNOG_OG1').is_in(ogs))
post_len = len(emapper)
print(f'Removed {pre_len - post_len:,} OGs (out of {pre_len:,}; {(pre_len - post_len)/pre_len:.2%})')

aa_sizes_by_og = defaultdict(list)
for og, unigene_len in emapper[['eggNOG_OG1', 'aa_size']].iter_rows():
    aa_sizes_by_og[og].append(unigene_len)

aa_sizes_by_og = {k:np.array(v) for k,v in aa_sizes_by_og.items()}
for v in aa_sizes_by_og.values():
    v.sort()


og_sizes = pl.DataFrame(
        [ # Multiply by 3 to get nucleotide length
            [k, v.mean() * 3, (v*3).std()]
            for k,v in aa_sizes_by_og.items()],
        orient='row',
        schema={
            'OG': pl.String,
            'mean': pl.Float32,
            'std': pl.Float32,
            }
        )


view_muts = [0, 100, 1000]
zscores_outs = []
all_zscores = []
all_esgs = []
esgs_fig,esgs_ax = plt.subplots()

for ifile, tag in jugspace['input_genomes']:

    emapper_out = jugspace['expanded'][tag]
    emapper_out = jug.value(emapper_out)
    emapper_out = emapper_out[['#query', 'eggNOG_OGs', 'original']]
    emapper_out = emapper_out.with_columns(
        eggNOG_OG1=emapper_out['eggNOG_OGs'].map_elements(
            extract_og,
            return_dtype=pl.String
            ))
    with gzip.open(jug.value(jugspace['gene_positions'][tag]), 'rb') as f:
        gene_positions = pl.read_csv(f, separator='\t')

    emapper_out = emapper_out.join(gene_positions, left_on='original', right_on='gene_id')[['nr_mutations', 'eggNOG_OG1', 'length']]
    emapper_out = emapper_out.filter(pl.col('eggNOG_OG1').is_in(ogs))
    all_muts = sorted(set(emapper_out['nr_mutations']))

    emapper_out = emapper_out.join(og_sizes, left_on='eggNOG_OG1', right_on='OG')

    zscores = emapper_out.with_columns(
        z=(pl.col('length')- pl.col('mean'))/pl.col('std')
        )[['nr_mutations', 'eggNOG_OG1', 'z']]


    esgs = []
    for mut in all_muts:
        sel = zscores.filter((pl.col('nr_mutations') == mut) & (pl.col('z') < -4)).select(pl.col('z').sum())
        esgs.append(sel.item())

    [ell] = esgs_ax.plot(all_muts, esgs, label=tag)
    esgs_ax.plot(all_muts, [esgs[0] for _ in all_muts], label=None, linestyle='--', color=ell.get_color())
    all_esgs.append(pl.DataFrame({'nr_mutations': all_muts, 'esg': esgs}).with_columns(tag=pl.lit(tag)))

    all_zscores.append(zscores.with_columns(tag=pl.lit(tag)))

    if PLOT_KO_LENGTHS_BOXPLOT:
        fig, ax = plt.subplots()
        sns.boxplot(data=zscores.filter(pl.col('nr_mutations').is_in(view_muts)),
                    x='nr_mutations',
                    y='z',
                    boxprops=dict(alpha=.4, facecolor='w'),
                    width=.2,
                    ax=ax,
                    )
        ax.set_title(f'z-scores for OG lengths ({tag})')
        fig.savefig(f'plots/emapper_ref/length_by_og_checkm2_zscore_{tag}.png')
        plt.close(fig)

    if PLOT_KO_LENGTHS_CUMMDIST:
        fig,ax = plt.subplots()
        X = np.linspace(-12,12,1000)
        for lim,lab in [(0, 'WT'), (100, '100'), (1000, '1k')]:
            sel = zscores.filter(pl.col('nr_mutations') == lim)
            ax.plot(X, [(sel['z'] < x).mean() for x in X], label=lab)
        ax.legend()
        ax.set_xlabel('OG length z-score')
        ax.set_ylabel('Fraction of genes (cumm)')
        ax.set_title(f'z-score distribution for OG lengths ({tag})')
        fig.tight_layout()
        fig.savefig(f'plots/emapper_ref/OG_zscores_{tag}.png')
        plt.close(fig)

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

all_esgs = pl.concat(all_esgs)
all_zscores = pl.concat(all_zscores)

esgs_ax.set_xlabel('Number of mutations')
esgs_ax.set_ylabel('Sum of z-scores < -4')
esgs_ax.legend()
sns.despine(esgs_fig, trim=True)
esgs_fig.tight_layout()
esgs_fig.savefig('plots/emapper_ref/esgs_m4_og.svg')
esgs_fig.savefig('plots/emapper_ref/esgs_m4_og.png')
