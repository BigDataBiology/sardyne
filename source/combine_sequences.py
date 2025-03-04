from glob import glob
from fasta import fasta_iter
from collections import defaultdict


def short_hash(seq):
    import hashlib
    return hashlib.md5(seq.encode()).hexdigest()[:10]

basedir = 'outputs/checkm2_103836.SAMN12024885_simulation/protein_files/'

sequence_presence = defaultdict(list)
sequences = set()

for f in glob(f'{basedir}/*.faa'):
    nr_mut = int(f.split('/')[-1].split('.')[0].split('_')[-1])
    for header, seq in fasta_iter(f):
        code = short_hash(seq)
        sequence_presence[code].append((nr_mut, header))
        sequences.add((code, seq))

with open(f'{basedir}/all_sequences.faa', 'wt') as f:
    for code, seq in sequences:
        f.write(f'>{code}\n{seq}\n')
