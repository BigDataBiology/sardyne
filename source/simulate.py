from jug import Task, TaskGenerator
from fasta import fasta_iter
import random
import gzip
import tempfile
import os

def prodigal_gene_sizes(input_file):
    '''Run prodigal and return the gene sizes (in nucleotides)'''
    import subprocess
    import numpy as np
    with tempfile.TemporaryDirectory() as tmpdirname:
        aa_out = f'{tmpdirname}/prodigal.faa'
        nt_out = f'{tmpdirname}/prodigal.fna'
        prodigal_out = f'{tmpdirname}/prodigal.txt'
        prodigal_err = f'{tmpdirname}/prodigal.err'
        zcat_p = subprocess.Popen(['zcat', input_file], stdout=subprocess.PIPE)
        with open(prodigal_err, 'w') as err:
            prodigal_p = subprocess.Popen(
                ['prodigal',
                     '-a', aa_out,
                     '-d', nt_out,
                     '-o', prodigal_out,
                 ],
                stderr=err,
                stdin=zcat_p.stdout
                )
        zcat_p.wait()
        prodigal_p.wait()
        gene_sizes = []
        for h, seq in fasta_iter(nt_out):
            gene_sizes.append(len(seq))
        return np.array(gene_sizes)

@TaskGenerator
def prodigal_gene_sizes_on_mut(input_file, nr_mutations):
    with tempfile.TemporaryDirectory() as tdir:
        create_mutated_file(f'{tdir}/mutated.fna.gz', input_file, nr_mutations)
        return prodigal_gene_sizes(f'{tdir}/mutated.fna.gz')

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
    headers = []
    for h, seq in fasta_iter(ifile):
        headers.append(h)
    assert len(headers) == 1
    return seq

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
            'conda', 'run', '-n', 'checkm2',
                'checkm2', 'predict',
                    '--threads', '8',
                     '--input', checkm2_idir,
                     '-x', 'fna.gz',
                     '--output-directory', odir,
                     ])
    return odir

INPUT_DATA = [
        ('ecoli_k12', '511145.SAMN02604091.fna.gz'),
        ('bacillus_subtilis', '1052585.SAMN02603352.fna.gz'),
        ('listeria_monocytogenes', '169963.SAMEA3138329.fna.gz'),
        ('campylobacter_jejuni', '192222.SAMEA1705929.fna.gz'),
        ('staphylococcus_aureus', '93061.SAMN02604235.fna.gz'),
        ]

nr_muts = list(range(0, 5_000, 25))

gene_sizes = {}
random_gene_sizes = {}
for tag, fname in INPUT_DATA:
    ifile = f'../data/genomes/{fname}'
    seq = read_seq(ifile)
    gene_sizes[tag] = []
    for i in nr_muts:
        gene_sizes[tag].append(prodigal_gene_sizes_on_mut(seq, i))
    run_checkm2(tag, seq, nr_muts)

    for method in ['uniform', 'markov2', 'markov4']:
        random_file = create_random_file(tag, seq, method)
        random_gene_sizes[tag, method] = Task(prodigal_gene_sizes, random_file)
