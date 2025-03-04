from glob import glob
import seaborn as sns
import polars as pl
from fasta import fasta_iter
from collections import defaultdict
from matplotlib import pyplot as plt


def short_hash(seq):
    import hashlib
    return hashlib.md5(seq.encode()).hexdigest()[:10]

basedir = 'outputs/checkm2_103836.SAMN12024885_simulation/protein_files/'

sequence_presence = defaultdict(list)
sequences = set()

sequences_per_nr_mut = defaultdict(set)

for f in glob(f'{basedir}/mutated*.faa'):
    nr_mut = int(f.split('/')[-1].split('.')[0].split('_')[-1])
    for header, seq in fasta_iter(f):
        code = short_hash(seq)
        sequence_presence[code].append((nr_mut, header))
        sequences.add((code, seq))
        sequences_per_nr_mut[nr_mut].add(code)

with open(f'{basedir}/all_sequences.faa', 'wt') as f:
    for code, seq in sequences:
        f.write(f'>{code}\n{seq}\n')

wt = sequences_per_nr_mut[0]
data = []
for m in sorted(sequences_per_nr_mut.keys()):
    cur = sequences_per_nr_mut[m]
    data.append((m, len(cur), len(cur - wt), len(wt - cur), len(cur & wt)))

data = pl.DataFrame(data, schema=['nr_mut', 'n', 'n_introduced', 'n_missing', 'n_common'], orient='row')
fig, ax = plt.subplots()
ax.plot('nr_mut', 'n_introduced', data=data, label='Introduced')
ax.plot('nr_mut', 'n_missing', data=data, label='Missing')
ax.plot('nr_mut', 'nr_mut', data=data, label=None, linestyle='--', color='black')
ax.set_xlabel('Number of mutations')
ax.legend()
fig.tight_layout()
sns.despine(fig, trim=True)
fig.savefig(f'plots/introduced_mutations.png')

