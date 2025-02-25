# Set up environment

We are using [pixi](https://pixi.sh/). Therefore, you should be able to activate the environment by running the following command:

```bash
pixi shell --environment sardyne
```

## Installing checkM2

*Note*: checkM2 requires Python &lt;3.9, so it is installed in a separate environment (`checkm2`). Before running it, you need to download the database.

```bash
pixi run --environment checkm2 checkm2 checkm2 database --download
```

This database is required for checkM2 to work properly.

