from jug import TaskGenerator

@TaskGenerator
def download_file(url, filename):
    import requests
    from os import makedirs
    makedirs('data', exist_ok=True)
    filename = f'data/{filename}'
    r = requests.get(url, stream=True)
    with open(filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
    return filename

@TaskGenerator
def dowload_from_progenomes(tax_id, sample_id):

    import requests
    import os
    URL = f'https://progenomes.embl.de/dumpSequence.cgi?p={tax_id}.{sample_id}&t=c&a={tax_id}'
    os.makedirs('genomes', exist_ok=True)
    r = requests.get(URL, stream=True)
    filename = f'genomes/{tax_id}.{sample_id}.fna.gz'
    with open(filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
    return filename

@TaskGenerator
def download_e_coli_k12():
    return dowload_from_progenomes.f('511145', 'SAMN02604091')

download_e_coli_k12()

# B. subtilis
dowload_from_progenomes('1052585', 'SAMN02603352')
# L. monocytogenes
dowload_from_progenomes('169963', 'SAMEA3138329')
# Campylobacter jejuni
dowload_from_progenomes('192222', 'SAMEA1705929')
# S. aureus
dowload_from_progenomes('93061', 'SAMN02604235')
# P. copri
dowload_from_progenomes('165179', 'SAMEA5853203')
# F. mortiferum
dowload_from_progenomes('469616', 'SAMN02463687')
# Streptococcus pneumoniae
dowload_from_progenomes('488222', 'SAMN02603444')

download_file('https://data.ace.uq.edu.au/public/misc_downloads/annotree/r83/uniref100.KO.faa', 'uniref100.KO.faa')

download_file('https://gmgc.embl.de/downloads/v1.0/subcatalogs/GMGC10.95nr.complete.faa.gz', 'GMGC10.95nr.complete.faa.gz')

download_file('https://swifter.embl.de/~fullam/spire/compiled/Coelho_2018_dog_spire_v1_MAGs.tar', 'Coelho_2018_dog_spire_v1_MAGs.tar')

download_file('https://progenomes.embl.de/data/repGenomes/progenomes3.contigs.representatives.fasta.bz2', 'progenomes3.contigs.representatives.fasta.bz2')
