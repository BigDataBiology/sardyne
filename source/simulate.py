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

INPUT_FILE = '../data/genomes/511145.SAMN02604091.fna.gz'
ecoli_k12 = read_seq(INPUT_FILE)
gene_sizes = []
nr_muts = list(range(0, 50_000, 50))
for i in nr_muts:
    mutated_file = create_mutated_file('ecoli_k12', ecoli_k12, i)
    gene_sizes.append(prodigal_gene_sizes(mutated_file))

random_gene_sizes = {}
for method in ['uniform', 'markov2', 'markov4']:
    random_file = create_random_file('ecoli_k12', ecoli_k12, method)
    random_gene_sizes[method] = prodigal_gene_sizes(random_file)


