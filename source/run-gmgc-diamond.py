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

chunks = split_fasta_file(GMGC_COMPLETE, 2_000_000)

ko_chunks = []
for ch in bvalue(chunks):
    ko_chunks.append(run_diamond(ch))

concatenate_files(ko_chunks, 'outputs/GMGC10.complete.diamond.out.gz')
