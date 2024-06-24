# Set up environment

```bash
conda create -n sardyne python=3.11 \
    jug polars matplotlib numpy seaborn requests \
    prodigal pyrodigal ipython
```

*Note*: checkM2 requires Python &lt;3.9, so we will install it in another environment.

## Installing checkM2

To install checkM2, use the following commands:

```bash
conda create -n checkm2 bioconda::checkm2
conda activate checkm2
checkm2 database --download
```

The `checkm2 database --download` command will download the latest version of the checkM2 database. This database is required for checkM2 to work properly.

