# Set up environment

We are using [pixi](https://pixi.sh/). Therefore, you should be able to activate the environment by running the following command:

```bash
pixi shell --environment sardyne
```

## Installing checkM2 & eggnogmapper

### checkM2

*Note*: checkM2 requires Python &lt;3.9, so it is installed in a separate environment (`checkm2`). Before running it, you need to download the database.

```bash
pixi run --environment checkm2 checkm2 checkm2 database --download
```

This database is required for checkM2 to work properly.

### eggnogmapper

Eggnogmapper does not require a separate environment, but it does require the database to be downloaded.

We will download it to `data/emapper-data`:

```bash
pixi run --environment sardyne

cd data
download_eggnog_data.py --data_dir emapper-data
```


