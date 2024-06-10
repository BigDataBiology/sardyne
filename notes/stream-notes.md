# Stream notes

## 2024-06-11

- Answer for why the previous results were wrong: CTRL-C does not always correctly stop checkM2
- [x] Generalize plotting code to multiple genomes

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

