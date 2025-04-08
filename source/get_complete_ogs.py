import polars as pl
from fasta import fasta_iter
from jug import Task

def extract_og(ogs):
    has_2 = None
    for c in ogs.split(','):
        og,taxon = c.split('@')
        if taxon == '1':
            return og
        elif taxon == '2':
            has_2 = og
    if has_2:
        return has_2
    if ogs: return ogs.split(',')[0].split('@')[0]


@Task
def complete_og_size_table():
    """
    Get the complete ogs from the GMGC10 dataset.
    """
    complete_sizes = []
    for h,seq in fasta_iter('../data/data/GMGC10.95nr.complete.faa.gz'):
        complete_sizes.append([h, len(seq)])

    complete_sizes = pl.DataFrame(complete_sizes, schema=['#query_name', 'aa_size'], orient='row')

    annotations = pl.scan_csv('../data/data/GMGC10.emapper2.annotations.tsv',
                    comment_prefix='# ',
                    separator='\t')
    annotations = annotations.select(['#query_name', 'eggNOG_OGs']).collect()
    annotations = annotations.with_columns(
            eggNOG_OG1=annotations['eggNOG_OGs'].map_elements(extract_og, return_dtype=str))
    annotations = annotations[['#query_name', 'eggNOG_OG1']]

    annotations = annotations.join(complete_sizes, on='#query_name', how='inner')
    oname = 'outputs/GMGC10.emapper2.annotations.complete.tsv'
    annotations.write_csv(oname,
                          separator='\t',
                          )
    return oname
