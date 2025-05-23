import polars as pl
from fasta import fasta_iter
from jug import TaskGenerator
from eggnog import extract_og

MIN_UNIGENES_PER_OG = 100


@TaskGenerator
def complete_og_size_table():
    """
    Get the complete ogs from the GMGC10 dataset.
    """
    import subprocess
    from os import path, unlink

    complete_sizes = []
    for h,seq in fasta_iter('../data/data/GMGC10.95nr.complete.faa.gz'):
        complete_sizes.append([h, len(seq)])

    complete_sizes = pl.DataFrame(complete_sizes, schema=['#query_name', 'aa_size'], orient='row')
    print("Complete sizes loaded")

    annotations = pl.scan_csv('../data/data/GMGC10.emapper2.annotations.tsv',
                    comment_prefix='# ',
                    separator='\t')
    annotations = annotations.with_columns(
            eggNOG_OG1=pl.col('eggNOG_OGs').map_elements(extract_og, return_dtype=str))
    annotations = annotations.select(['#query_name', 'eggNOG_OG1'])
    print("Annotations loaded (lazy)")

    annotations = annotations.join(complete_sizes.lazy(), on='#query_name', how='inner')
    print("Annotations joined with complete sizes")

    oname = 'outputs/GMGC10.emapper2.annotations.complete.tsv'
    annotations.sink_csv(oname,
                           separator='\t',
                           )
    if path.exists(oname + '.gz'):
        unlink(oname + '.gz')
    subprocess.check_call(['gzip', oname])
    return oname + '.gz'


@TaskGenerator
def baseline_og_sizes(gmgc10_complete, min_unigenes_per_og):
    import numpy as np
    import gzip
    from collections import defaultdict
    with gzip.open(gmgc10_complete, 'rb') as f:
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
    oname = 'outputs/GMGC10.emapper2.annotations.complete.ogsizes.tsv'
    og_sizes.write_csv(oname, separator='\t')
    return oname


baseline_og_sizes(complete_og_size_table(), MIN_UNIGENES_PER_OG)
