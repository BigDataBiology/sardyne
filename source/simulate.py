from jug import Task, TaskGenerator, bvalue
from fasta import fasta_iter
import random
import gzip
import tempfile
import re
import os

NR_CHECKM2_THREADS = 8

def select_every(g_id: str, n: int):
    import hashlib
    h = hashlib.sha256()
    h.update(b'sardyne')
    h.update(g_id.encode('ascii'))
    h = int(h.hexdigest(), 16)
    return h % n == 0


@TaskGenerator
def expand_progenomes3(n):
    os.makedirs('../data/data/progenomes3.expanded', exist_ok=True)

    prev = 'x'
    seen = set()
    for h, seq in fasta_iter('../data/data/progenomes3.contigs.representatives.fasta.bz2'):
        tokens = h.split('.')
        assert len(tokens) == 3
        g_id = '.'.join(tokens[:2])
        if g_id != prev:
            prev = g_id
            if select_every(g_id, n):
                oname = f'../data/data/progenomes3.expanded/{g_id}.fna.gz'
                out = gzip.open(oname, 'wt')
                assert oname not in seen
                seen.add(oname)
            else:
                out = None
        if out is not None:
            out.write(f'>{h}\n{seq}\n')
    return sorted(seen)


def mutate1(seq : list[str]) -> None:
    '''Mutate a sequence in place using a simple model'''
    pos = random.randint(0, len(seq) - 1)
    op_r = random.random()
    if op_r < .4:
        # Insertion
        nc = random.choice('ACGT')
        seq[pos] += nc
    elif op_r < .8:
        # Deletion
        while seq[pos] == '':
            pos = random.randint(0, len(seq) - 1)
        seq[pos] = seq[pos][1:]
    else:
        # Substitution
        nc = seq[pos]
        while nc == seq[pos]:
            nc = random.choice('ACGT')
        seq[pos] = nc

def mutate_multi(seq, n=1):
    '''Mutate a sequence n times'''
    random.seed(seq[:1024])
    seq = list(seq)
    for i in range(n):
        mutate1(seq)
    return ''.join(seq)

def create_mutated_file(oname, seq, n):
    '''Create a mutated file with n mutations'''
    seq = mutate_multi(seq, n)
    with gzip.open(oname, 'wt', compresslevel=0) as f:
        f.write(f'>mutated_{n}\n{seq}\n')
    return oname

@TaskGenerator
def read_seq(ifile):
    '''Read a sequence from a FASTA file, checking that there is only one sequence'''
    seqs = []
    for _, seq in fasta_iter(ifile):
        seqs.append(seq)
    padding = ''.join(['N' for _ in range(300)])
    return padding.join(seqs)

def random_dna_same_len_as(seq):
    return ''.join(random.choices('ACGT', k=len(seq)))

def random_dna_same_len_as_markov_chain(seq, mc_len=2):
    '''Generate a random DNA sequence with the same length and Markov chain as seq'''
    import numpy as np
    from collections import Counter
    import random
    counts_mc_len = Counter(seq[i:i+mc_len] for i in range(len(seq)-mc_len))
    counts_mc_lenp1 = Counter(seq[i:i+mc_len+1] for i in range(len(seq)-mc_len-1))
    probs = {}
    for k,v in counts_mc_len.items():
        k = tuple(k)
        probs[k] = []
        for n in 'ATCG':
            probs[k].append((counts_mc_lenp1.get(''.join(k)+n, 0) + 1) / (v + 4))
    new_seq = []
    for _ in range(mc_len):
        new_seq.append(random.choice('ATCG'))
    for i in range(len(seq)-mc_len):
        new_seq.append(random.choices('ATCG', probs[tuple(seq[-mc_len:])])[0])
    return ''.join(new_seq)

@TaskGenerator
def create_random_file(tag, seq, method):
    '''Create a random file with the same length as seq'''
    os.makedirs('outputs/simulations', exist_ok=True)
    ofile = f'outputs/simulations/{tag}_random_{method}.fna.gz'
    if method == 'markov2':
        seq = random_dna_same_len_as_markov_chain(seq, 2)
    elif method == 'markov4':
        seq = random_dna_same_len_as_markov_chain(seq, 4)
    elif method == 'uniform':
        seq = random_dna_same_len_as(seq)
    else:
        raise NotImplementedError(f'Unknown method {method}')

    seq = random_dna_same_len_as(seq)
    with gzip.open(ofile, 'wt') as f:
        f.write(f'>{tag}_random\n{seq}\n')
    return ofile

@TaskGenerator
def run_checkm2(tag, seq, nr_muts):
    import subprocess
    import shutil
    import tempfile
    import os
    import pathlib
    with tempfile.TemporaryDirectory() as tdir:
        checkm2_idir = f'{tdir}/checkm2_inputs'
        os.makedirs(checkm2_idir)
        for i in nr_muts:
            create_mutated_file(f'{checkm2_idir}/mutated_{i:06}.fna.gz', seq, i)
        odir = f'outputs/checkm2_{tag}_simulation'
        shutil.rmtree(odir, ignore_errors=True)
        subprocess.check_call([
            'pixi', 'run', '--environment', 'checkm2',
                'checkm2', 'predict',
                    '--threads', str(NR_CHECKM2_THREADS),
                     '--input', checkm2_idir,
                     '-x', 'fna.gz',
                     '--output-directory', odir,
                     ])
    return odir


PRODIGAL_FASTA_HEADER_PAT = r'^>(\S+) # (\d+) # (\d+) # (-?)1 # '
def get_gene_positions(f: str):
    '''Extract gene positions from a Prodigal FASTA file'''
    import polars as pl
    data = []
    for line in open(f):
        if match := re.match(PRODIGAL_FASTA_HEADER_PAT, line):
            gene_id, start, end, is_reverse = match.groups()
            start = int(start)
            end = int(end) + 1
            is_reverse = is_reverse == '-'
            length = abs(end - start)
            data.append((gene_id, start, end, length, is_reverse))
    return pl.DataFrame(data, schema={
                    'gene_id': pl.String,
                    'start': pl.Int32,
                    'end': pl.Int32,
                    'length': pl.Int32,
                    'is_reverse': pl.Boolean},
                    orient='row')

@TaskGenerator
def get_all_gene_positions(c2_odir, nr_muts):
    from glob import glob
    import polars as pl
    partials = []

    faas = sorted(glob(f'{c2_odir}/protein_files/*.faa'))
    assert len(faas) == len(nr_muts)
    for f in faas:
        match = re.match(r'.*mutated_(\d+).fna.faa', f)
        nr_mut = int(match.group(1))
        partials.append(
            get_gene_positions(f).with_columns(
            nr_mutations=pl.lit(nr_mut)
            ))
    positions = pl.concat(partials)
    oname = f'{c2_odir}/all_gene_positions.tsv.gz'
    with gzip.open(f'{c2_odir}/all_gene_positions.tsv.gz', 'wb') as f:
        positions.write_csv(f, separator='\t')

    # Not 100% safe as we can get to a state where there is an
    # interruption in the cleanup loop, leaving the process in a bad
    # state
    for f in faas:
        os.unlink(f)
    return oname


INPUT_DATA = [
        ('ecoli_k12', '511145.SAMN02604091.fna.gz'),
        ('bacillus_subtilis', '1052585.SAMN02603352.fna.gz'),
        ('listeria_monocytogenes', '169963.SAMEA3138329.fna.gz'),
        ('campylobacter_jejuni', '192222.SAMEA1705929.fna.gz'),
        ('staphylococcus_aureus', '93061.SAMN02604235.fna.gz'),
        ('prevotella_copri', '165179.SAMEA5853203.fna.gz'),
        ('fusobacterium_mortiferum', '469616.SAMN02463687.fna.gz'),
        ('streptococcus_pneumoniae', '488222.SAMN02603444.fna.gz'),
        ]

nr_muts = list(range(0, 5_000, 25))

input_genomes = bvalue(expand_progenomes3(800))
for ifile in input_genomes:
    tag = ifile.split('/')[-1].removesuffix('.fna.gz')
    seq = read_seq(ifile)
    c2_odir = run_checkm2(tag, seq, nr_muts)
    get_all_gene_positions(c2_odir, nr_muts)
