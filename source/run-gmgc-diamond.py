from jug import TaskGenerator, bvalue
NR_THREADS = 8
GMGC_COMPLETE = '../data/data/GMGC10.95nr.complete.faa.gz'

@TaskGenerator
def run_diamond(ifile):
    import subprocess
    from pathlib import Path
    from os import unlink
    ofile = ifile.replace('.faa.gz', '.tsv')
    subprocess.check_call(
        ['conda', 'run', '-n', 'checkm2',
             'diamond', 'blastp',
             '--outfmt', '6',
             '--max-target-seqs', '1',
             '--query', ifile,
             '-o', ofile,
             '--threads', str(NR_THREADS),
             '--db', (Path.home() / 'databases/CheckM2_database/uniref100.KO.1.dmnd'),
             '--query-cover', '80',
             '--subject-cover', '80',
             '--id', '30',
             '--evalue', '1e-05',
             '--block-size', '2.0'])
    unlink(ifile)
    return ofile

@TaskGenerator
def split_fasta_file(ifile, nr_seqs_per_chunk):
    import gzip
    from os import makedirs
    import fasta
    makedirs('outputs/gmgc_diamond', exist_ok=True)

    ofiles = []
    cur = nr_seqs_per_chunk # to trigger a new file at start
    ix = 0
    ofile = None
    for h,seq in fasta.fasta_iter(ifile):
        if cur == nr_seqs_per_chunk:
            if ofile is not None:
                ofile.close()
            cur = 0
            ofiles.append(f'outputs/gmgc_diamond/chunk_{ix:02}.faa.gz')
            ofile = gzip.open(ofiles[-1], 'wt', compresslevel=1)
            ix += 1
        cur += 1
        ofile.write(f'>{h}\n{seq}\n')
    ofile.close()
    return ofiles

@TaskGenerator
def concatenate_files(ofiles, ofile):
    import gzip
    from os import unlink
    with gzip.open(ofile, 'wb') as of:
        for f in ofiles:
            with open(f, 'rb') as f:
                while ch := f.read(1024*1024):
                    of.write(ch)
    for f in ofiles:
        unlink(f)
    return ofile

@TaskGenerator
def annotate_with_orf_lengths(ifile, fareference, ofile):
    import gzip
    import polars as pl
    from fasta import fasta_iter
    with gzip.open(ifile, 'rb') as f:
        gmgc = pl.read_csv(f,
                           separator='\t',
                           has_header=False,
                           columns=[0, 1, 2, 3, 4, 6, 7, 8, 9, ],
                           new_columns=[
                               'query', # 0
                               'subject', # 1
                               'identity', # 2
                               'alignment_length', # 3
                               'mismatches', # 4
                               # 'gap_opens', # 5
                               'q_start', # 6
                               'q_end', # 7
                               's_start', # 8
                               's_end', # 9
                               # 'evalue', # 10
                               # 'bit_score', # 11
                               ])

    interesting_unigenes = set(gmgc['query'].to_list())
    unigene_to_len = {}
    for h,seq in fasta_iter(fareference):
        if h in interesting_unigenes:
            unigene_to_len[h] = len(seq)

    gmgc = gmgc.with_columns([
        gmgc['query'].map_elements(unigene_to_len.get, return_dtype=int).alias('gmgc_unigene_len_aa'),
        ])

    with gzip.open(ofile, 'wt') as of:
        gmgc.write_csv(of, separator='\t')
    return ofile

chunks = split_fasta_file(GMGC_COMPLETE, 2_000_000)

ko_chunks = []
for ch in bvalue(chunks):
    ko_chunks.append(run_diamond(ch))

concat_outs = concatenate_files(ko_chunks, 'outputs/GMGC10.complete.diamond.out.gz')
annotate_with_orf_lengths(concat_outs,
                          GMGC_COMPLETE,
                          'outputs/GMGC10.complete.annotated.tsv.gz')
