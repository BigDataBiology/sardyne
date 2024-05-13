from jug import TaskGenerator

@TaskGenerator
def download_e_coli_k12():
    import requests
    import os
    URL = 'https://progenomes.embl.de/dumpSequence.cgi?p=511145.SAMN02604091&t=c&a=511145'
    r = requests.get(URL, stream=True)
    filename = 'genomes/511145.SAMN02604091.fna.gz'
    with open(filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
    return filename

download_e_coli_k12()
