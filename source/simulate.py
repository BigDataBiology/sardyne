from jug import TaskGenerator
from fasta import fasta_iter
import random
import gzip


@TaskGenerator
def prodigal_gene_sizes(input_file):
    '''Run prodigal and return the gene sizes (in nucleotides)'''
    import subprocess
    import tempfile
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

def mutate1(seq):
    pos = random.randint(0, len(seq) - 1)
    op_r = random.random()
    if op_r < .4:
        # Insertion
        nc = random.choice('ACGT')
        seq = seq[:pos] + nc + seq[pos:]
    elif op_r < .8:
        # Deletion
        seq = seq[:pos] + seq[pos+1:]
    else:
        # Substitution
        nc = seq[pos]
        while nc == seq[pos]:
            nc = random.choice('ACGT')
        seq = seq[:pos] + nc  + seq[pos+1:]
    return seq

def mutate_multi(seq, n=1):
    '''Mutate a sequence n times'''
    random.seed(seq[:1024])
    for i in range(n):
        seq = mutate1(seq)
    return seq

@TaskGenerator
def create_mutated_file(tag, seq, n):
    '''Create a mutated file with n mutations'''
    ofile = f'outputs/simulations/{tag}_mutated_{n}.fna.gz'
    seq = mutate_multi(seq, n)
    with gzip.open(ofile, 'wt') as f:
        f.write(f'>{tag}_mutated_{n}\n{seq}\n')
    return ofile

@TaskGenerator
def read_seq(ifile):
    '''Read a sequence from a FASTA file, checking that there is only one sequence'''
    headers = []
    for h, seq in fasta_iter(ifile):
        headers.append(h)
    assert len(headers) == 1
    return seq

INPUT_FILE = '../data/genomes/511145.SAMN02604091.fna.gz'
ecoli_k12 = read_seq(INPUT_FILE)
gene_sizes = []
nr_muts = list(range(0, 50_000, 50))
for i in nr_muts:
    mutated_file = create_mutated_file('ecoli_k12', ecoli_k12, i)
    gene_sizes.append(prodigal_gene_sizes(mutated_file))


