# Stream notes

## Next streams

- Use the dog gut MAGs from Coelho et al., 2018 (eventually expand to other datasets)
- Attempt to find KOs with very stable gene sizes
- Look at all KO hits (not just those that are unique to a genome)
- Try filtering DIAMOND hits by position (i.e., only keep those that align near the start of the genes)
- Explore variations of the ESGS metric (for example, what is the best threshold for the z-score?)

## 2025-03-18

- Moved code to de-replicate sequences into `simulate.py` (and into the `main` branch)

## 2025-03-04

Live: https://www.youtube.com/watch?v=YDNmvQpQ2DI

- Added eggnogmapper and downloaded its database
- Explored how gene sequences change with (simulated) mutations (in `wip` branch)


## 2025-02-25

Recap: https://youtu.be/1GV5VpvNkTI
Live: https://www.youtube.com/watch?v=3EibnfuPUMc

- Updated to pixi (from conda)
- Cleaned plotting code to avoid creating too many plots
- Browsed some outliers
- Added the special case microbes (e.coli et al) back into the analysis

## 2025-02-18

Recap: https://youtu.be/8E0DTEPSUhw
Live: https://youtube.com/live/cLDxd6EyLPo

- Cleaned code to remove redundancy
- Save fewer intermediate files to enable scaling up
- Use a selection of representative genomes from ProGenomes3 (1%)

## 2025-02-11

Recap: https://youtu.be/OE6BSEDDqd0
Live: https://www.youtube.com/live/C5zeyqYQGbg

- Cleaned up the code from last week & merged to `main`
- Expanded the dataset to a few more hand-picked species, but also downloaded the full ProGenomes3 representative set
- Downloaded dog gut MAGs from [Coelho et al., 2018](https://dx.doi.org/10.1186/s40168-018-0450-3)


## 2025-02-04

Recap: https://youtu.be/iXQRjUq4kAE
Live: https://www.youtube.com/live/rQp-OkmpHxY

- Proposed metric: sum of z-scores for gene sizes per KO below a certain threshold (e.g., -4). Proposed name: excess of short genes score (ESGS)
- Plotted ESGS for a few example species

## 2025-01-28

- Since _E. coli_ is part of the database, we recover the exact genes
- K02436 seems to have a relatively stable gene size in the GMGCv1 complete gene set (227-304 AA, with the mean 254 and std. dev of 8.2), but the E. coli hit is much shorter (148 AA)
- K01232 is less clearly an outlier, but the GMGCv1 smallest hit **is the E. coli** hit at 212 (otherwise, mean is 443 and std. dev. 18.4)

## 2024-10-01

- Cleanup the `length_by_ko.py` script: (1) use all ORF/KOs and (2) refactor and optimize the code

## 2024-09-16

- Precompute GMGCv1 KO sizes to avoid rerunning every time (integrated with run-gmgc-diamond Jug script)
- Cleanup `length_by_ko.py` script

## 2024-09-09

- Recap: https://youtu.be/7HE_au3jYSY
- Full stream: https://www.youtube.com/live/VwvZ5DEXV-E

- Used [GMGCv1](https://gmgc.embl.de/) as database for KO sizes
- Concatenated GMGCv1 DIAMOND outputs into a single file

## 2024-09-03

- Recap: https://youtu.be/DDfLfcmnpVI
- Full stream: https://www.youtube.com/live/LLeJDXkNFBY

- Split [GMGCv1](https://gmgc.embl.de/) into chunks of 2m sequences to make the KO mapping more robust
- Checked zscores per KO for gene sizes (instead of percentiles): appears much more stable across different species.

## 2024-08-13

- Recap: https://youtu.be/a0P3La0w91g
- Full stream: https://www.youtube.com/live/h_BpiVLNy7M

- Mapped GMGCv1 to KOs using diamond (same version as checkM2)
- Explored data. Looked at the distribution of increase in gene sizes vs. decrease. It is not symmetric as it is easier to introduce a stop codon than to read through random DNA without hitting a stop codon sooner or later.

## 2024-08-06

- Recap: https://youtu.be/d2htGMZOVBY
- Full stream: https://youtube.com/live/anSEVzkygIU

- Downloaded GMGCv1 database (complete genes only, amino acid sequences)
- Checked percentiles of gene sizes by KO using checkM2 reference itself

## 2024-06-25

- Recap: https://youtu.be/LmlGyFV7gXk
- Full stream: https://www.youtube.com/watch?v=c9vHbBqc1R0

- Used Uniref100.KO database to get estimates of sizes by KO (work committed to `wip` branch)
- Installed seaborn into sardyne conda environment

## 2024-06-18

- Recap video: https://youtu.be/21m5qY9ovKw
- Full stream: https://www.youtube.com/live/ivIOWYNPmMg

- Explored length of ORFs by KO in the E. coli simulations [commited to `wip` branch]
- Changed the number of mutations in simulations to max 5k by steps of 25

## 2024-06-11

- Recap video: https://youtu.be/hYwsaFxq5dc
- Full stream: https://www.youtube.com/live/uSTN6O44NCg

- Answer for why the previous results were wrong: CTRL-C does not always correctly stop checkM2
- [x] Generalize plotting code to multiple genomes
- [x] Analysed checkM2 results by looking at missing/false positive KOs in mutated genomes

## 2024-06-04

- Recap video: https://youtu.be/nQCxzccHhnY
- Full stream: https://www.youtube.com/live/aaU-NqOID2w

Plans for stream

- [x] Plot checkM2 results on simulated data
- [x] Explore a different genome (high-quality isolate)
- check gene sizes on lower-quality isolates from ProGenomes


## 2024-05-28

- Recap video: https://youtu.be/RjP3IaKr0Mo (extra on [Jug](https://jug.rtfd.io/): https://youtu.be/NYiEE6ok9Ds
- Full stream: https://youtube.com/live/WONbVZtHG64

- [x] run checkM2 (created code to run checkM2 on the simulated data from jug)
- [x] create a completely noise genome (random sequence with same length as the E. coli genome)

## 2024-05-21

- Recap video: https://youtu.be/k_aaJn_FP1g
- Full stream: https://www.youtube.com/live/WBGHgXiZFM4

- Simulation of mutations (with a very simple, even simplistic model)
- Measure average gene size (predicted by prodigal)
- Used checkM2 as well

Next steps: use checkM2 on the simulated data to get its estimates of the genomic quality.

## 2024-05-14

- Recap video: https://youtu.be/N-O3n63-fcY
- Full stream (with audio issues): https://www.youtube.com/live/D4HINfXUGgA

1. Presenting the problem
2. Set up repository/environment

### References

Many purported pseudogenes in bacterial genomes are bona fide genes
Nicholas P. Cooley & Erik S. Wright
https://link.springer.com/article/10.1186/s12864-024-10137-0


CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning
Alex Chklovski, Donovan H. Parks, Ben J. Woodcroft & Gene W. Tyson
https://www.nature.com/articles/s41592-023-01940-w

## 2024-05-06

Pre-stream notes (shared on the mailing-list):

I will go into more detail next week on the stream, but the project for season 2 of EOS is to develop a tool for estimating genome quality using pseudo-genes.

Some background:

When using Nanopore (ONT, for Oxford Nanopore Technology), we can get genomes (MAGs, using, for example, SemiBin2) that are very complete, but have a high rate of SNPs/indels. This results in a lot of apparent pseudo-genes, but many are probably artefacts. A recent paper also argued this:

>  Many purported pseudogenes in bacterial genomes are bona fide genes by Cooley & Wright
>  http://dx.doi.org/10.1186/s12864-024-10137-0

When we use tools like checkM2, this gets picked up indirectly as missing genes and a loss of completeness, but the idea is to explore more directly whether we can look for an over-abundance of apparent pseudo-genes (pseudo pseudo genes?) to estimate quality.

The basic outline is to develop a tool to do something like the following:

1. take the MAG and call genes using prodigal (or pyrodigal)
2. annotate these genes using eggnog-mapper to assign them to OGs
3. check if the genes are smaller than expected
4. convert the outputs into a quality score
5. This is all still very hand-waved, but this is the basic outline.

We might start first by looking into the problem by simulating data and seeing how genes get shorter/quality drops as we introduce more errors.

