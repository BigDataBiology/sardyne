from jug import TaskGenerator, bvalue
from fasta import fasta_iter
import random
import gzip
import re
import os
from os import path
from eggnog import extract_og

NR_CHECKM2_THREADS = 8
NR_EMAPPER_THREADS = 8


def select_every(g_id: str, n: int):
    import hashlib
    h = hashlib.sha256()
    h.update(b'sardyne')
    h.update(g_id.encode('ascii'))
    h = int(h.hexdigest(), 16)
    return h % n == 0


@TaskGenerator
def expand_progenomes3(n, allow_list):
    os.makedirs('../data/data/progenomes3.expanded', exist_ok=True)

    allow_list = {g_id:tag for tag, g_id in allow_list}

    prev = 'x'
    seen = set()
    tags = {}
    for h, seq in fasta_iter('../data/data/progenomes3.contigs.representatives.fasta.bz2'):
        tokens = h.split('.')
        assert len(tokens) == 3
        g_id = '.'.join(tokens[:2])
        if g_id != prev:
            prev = g_id
            if select_every(g_id, n) or g_id in allow_list:
                oname = f'../data/data/progenomes3.expanded/{g_id}.fna.gz'
                if path.exists(oname):
                    out = None
                else:
                    out = gzip.open(oname, 'wt')
                assert oname not in seen
                seen.add(oname)
                if g_id in allow_list:
                    tags[oname] = allow_list[g_id]
                else:
                    tags[oname] = g_id
            else:
                out = None
        if out is not None:
            out.write(f'>{h}\n{seq}\n')
    seen = sorted(seen)
    return [(oname, tags[oname]) for oname in seen]


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
    from collections import Counter
    import random
    counts_mc_len = Counter(seq[i:i+mc_len] for i in range(len(seq)-mc_len))
    counts_mc_lenp1 = Counter(seq[i:i+mc_len+1] for i in range(len(seq)-mc_len-1))
    probs = {}
    for k, v in counts_mc_len.items():
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


PRODIGAL_FASTA_HEADER_PAT = r'(\S+) # (\d+) # (\d+) # (-?)1 # '
def get_gene_positions(headers):
    '''Extract gene positions from a list of headers from Prodigal'''
    import polars as pl
    data = []
    for line in headers:
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
def get_all_gene_positions(c2_odir, sequences_per_nr_mut):
    import polars as pl
    partials = []

    for nr_mut, c_headers in sequences_per_nr_mut.items():
        headers = [h for _, h in c_headers]
        partials.append(
            get_gene_positions(headers).with_columns(
                nr_mutations=pl.lit(nr_mut)
            ))
    positions = pl.concat(partials)
    oname = f'{c2_odir}/all_gene_positions.tsv.gz'
    with gzip.open(oname, 'wb') as f:
        positions.write_csv(f, separator='\t')

    return oname


def short_hash(seq):
    import hashlib
    return hashlib.md5(seq.encode()).hexdigest()[:24]


@TaskGenerator
def collate_protein_sequences(basedir, nr_muts):
    from collections import defaultdict
    from glob import glob

    sequences = set()
    sequences_per_nr_mut = defaultdict(set)

    faas = sorted(glob(f'{basedir}/protein_files/mutated*.faa'))
    assert len(faas) == len(nr_muts)
    for f in faas:
        nr_mut = int(f.split('/')[-1].split('.')[0].split('_')[-1])
        for header, seq in fasta_iter(f, full_header=True):
            code = short_hash(seq)
            sequences.add((code, seq))
            sequences_per_nr_mut[nr_mut].add((code, header))

    assert len(set(code for code, _ in sequences)) == len(sequences)

    oname = f'{basedir}/all_proteins.faa'
    with open(oname, 'wt') as f:
        for code, seq in sorted(sequences):
            f.write(f'>{code}\n{seq}\n')
    return oname, dict(sequences_per_nr_mut)


@TaskGenerator
def run_emapper(basedir, protein_file):
    import subprocess
    from pathlib import Path
    output = f'{basedir}/emapper'

    subprocess.check_call([
        'emapper.py',
            '-i', protein_file,
            '--output', output,
            '--cpu', str(NR_EMAPPER_THREADS),
            '--data_dir', Path.home() / 'databases/emapper_db',
            '--itype', 'proteins',
            '--tax_scope', 'auto',
            ])
    return output


@TaskGenerator
def cleanup_faas(basedir, nr_muts, run_after):
    from glob import glob
    faas = sorted(glob(f'{basedir}/protein_files/mutated*.faa'))
    assert len(faas) == len(nr_muts)
    for f in faas:
        os.unlink(f)


@TaskGenerator
def expand_emapper(emapper_outprefix, seqmaps):
    import polars as pl
    emapper_out = f'{emapper_outprefix}.emapper.annotations'
    emapper_out = pl.read_csv(emapper_out,
                              skip_lines=4,
                              separator='\t',
                              comment_prefix='##')
    hash2original = []
    for gs in seqmaps.values():
        for seqhash, header in gs:
            header, _ = header.split(' ', 1)
            hash2original.append((seqhash, header))

    hash2original = pl.DataFrame(hash2original,
                                 schema=['hash', 'original'],
                                 orient='row')
    return emapper_out.join(hash2original, left_on='#query', right_on='hash')


@TaskGenerator
def plot_simulation_results(expanded, tag, use_og1s=False):
    import polars as pl
    from matplotlib import pyplot as plt
    import seaborn as sns
    print(f'plotting {tag}')

    os.makedirs('plots/introduced_mutations', exist_ok=True)
    expanded = expanded.with_columns(
            nr_mutations=expanded['original']
                        .str.split('_')
                        .list[1]
                        .cast(pl.Int32)
                        .alias('nr_mutations'))

    if use_og1s:
        expanded = expanded.with_columns(
                    eggNOG_OG1=pl.col('eggNOG_OGs')
                            .map_elements(extract_og,
                            return_dtype=str))
        col = 'eggNOG_OG1'
    else:
        col = '#query'

    expanded = expanded[[col, 'nr_mutations']]

    wt = set(expanded.filter(expanded['nr_mutations'] == 0)[col])

    data = []
    for row in expanded.group_by('nr_mutations').all().iter_rows():
        m, seqs = row
        cur = set(seqs)
        data.append((m, len(cur), len(cur - wt), len(wt - cur), len(cur & wt)))

    data = pl.DataFrame(data, schema=['nr_mut', 'n', 'n_introduced', 'n_missing', 'n_common'], orient='row')
    data = data.sort('nr_mut')
    fig, ax = plt.subplots()
    ax.plot('nr_mut', 'n_introduced', data=data, label='Introduced')
    ax.plot('nr_mut', 'n_missing', data=data, label='Missing')
    if not use_og1s:
        ax.plot('nr_mut', 'nr_mut', data=data, label=None, linestyle='--', color='black')
    ax.set_xlabel('Number of mutations')
    ax.legend()
    fig.tight_layout()
    sns.despine(fig, trim=True)
    suffix = '_og1s' if use_og1s else ''
    oname = f'plots/introduced_mutations/{tag}{suffix}.png'
    fig.savefig(oname, dpi=300)
    return oname


@TaskGenerator
def select_columns(df, columns):
    '''Select columns from a DataFrame'''
    return df.select(columns)


@TaskGenerator
def load_og_sizes():
    '''Load OG sizes from a file'''
    import polars as pl
    return pl.read_csv('outputs/GMGC10.emapper2.annotations.complete.ogsizes.tsv', separator='\t')


@TaskGenerator
def zscore_by_og(emapper_out, gene_positions, og_sizes):
    '''Z-score the lengths of the genes by OG'''
    import polars as pl
    ogs = set(og_sizes['OG'])

    emapper_out = emapper_out[['#query', 'eggNOG_OGs', 'original']]
    emapper_out = emapper_out.with_columns(
        eggNOG_OG1=emapper_out['eggNOG_OGs'].map_elements(
            extract_og,
            return_dtype=pl.String
            ))
    with gzip.open(gene_positions, 'rb') as f:
        gene_positions = pl.read_csv(f, separator='\t')

    emapper_out = emapper_out.join(gene_positions, left_on='original', right_on='gene_id')[['nr_mutations', 'eggNOG_OG1', 'length']]
    emapper_out = emapper_out.filter(pl.col('eggNOG_OG1').is_in(ogs))

    emapper_out = emapper_out.join(og_sizes, left_on='eggNOG_OG1', right_on='OG')

    zscores = emapper_out.with_columns(
        z=(pl.col('length')- pl.col('mean'))/pl.col('std')
        )[['nr_mutations', 'eggNOG_OG1', 'z']]
    return zscores


SPECIAL_MICROBES = [
        ('ecoli_k12', '511145.SAMN02604091'),
        ('bacillus_subtilis', '1052585.SAMN02603352'),
        ('listeria_monocytogenes', '169963.SAMEA3138329'),
        ('campylobacter_jejuni', '192222.SAMEA1705929'),
        ('staphylococcus_aureus', '93061.SAMN02604235'),
        ('prevotella_copri', '165179.SAMEA5853203'),
        ('fusobacterium_mortiferum', '469616.SAMN02463687'),
        ('streptococcus_pneumoniae', '488222.SAMN02603444'),
        ]

nr_muts = list(range(0, 5_000, 25))

input_genomes = bvalue(expand_progenomes3(200, SPECIAL_MICROBES))
expanded = {}
gene_positions = {}
og_sizes = load_og_sizes()
zscores = {}
for ifile, tag in input_genomes:
    seq = read_seq(ifile)
    c2_odir = run_checkm2(tag, seq, nr_muts)
    oname_seqmaps = collate_protein_sequences(c2_odir, nr_muts)
    cleanup_faas(c2_odir, nr_muts, run_after=oname_seqmaps)

    emapper_outprefix = run_emapper(c2_odir, oname_seqmaps[0])
    expanded[tag] = expand_emapper(emapper_outprefix, oname_seqmaps[1])
    expanded[tag] = select_columns(expanded[tag], ['#query', 'eggNOG_OGs', 'original'])
    plot_simulation_results(expanded[tag], tag)
    plot_simulation_results(expanded[tag], tag, use_og1s=True)
    gene_positions[tag] = get_all_gene_positions(c2_odir, oname_seqmaps[1])
    zscores[tag] = zscore_by_og(expanded[tag], gene_positions[tag], og_sizes)

