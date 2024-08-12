from jug import TaskGenerator
NR_THREADS = 8

@TaskGenerator
def run_gmgc_diamond():
    import subprocess
    from pathlib import Path

    ofile = 'outputs/gmgc_diamond.tsv'

    subprocess.check_call(
        ['conda', 'run', '-n', 'checkm2',
             'diamond', 'blastp',
             '--outfmt', '6',
             '--max-target-seqs', '1',
             '--query', '../data/data/GMGC10.95nr.complete.faa.gz',
             '-o', ofile,
             '--threads', str(NR_THREADS),
             '--db', (Path.home() / 'databases/CheckM2_database/uniref100.KO.1.dmnd'),
             '--query-cover', '80',
             '--subject-cover', '80',
             '--id', '30',
             '--evalue', '1e-05',
             '--block-size', '2.0'])
    return ofile

run_gmgc_diamond()
