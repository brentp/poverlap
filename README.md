Poverlap
========
Simple, flexible, parallized significance testing of a pair of BED files.

BEDTools offers the basic machinery to perform repeated randomization and
overlap testing. 

This script automatically parallelizes that process. It has these additional
features:

 1) allows shuffling each interval to a location within N bases of its current
    location. This is a very basic start toward something like Bickels work for
    ENCODE that can preserve local structure

 2) given the `a` and `b` BED files specified in order, the default is to
    shuffle only the `b` BED file. However, both files may be shuffled with
    `shuffle_both` argument to poverlap.

 3) This implements the sampling schema defined in Haminen et al in BMC
    Bioinformatics 9: 336. Whereby we give a list of possible locations,
    e.g. the locations of a set of transcriptions factors and we wish to
    see how 2 of them e.g. CTCF and Pol2 are related. To do this, we fix
    the locations of CTCF (from within the BED file) and we randomize the
    location of Pol2 to any location occupied by a TF that is not CTCF.
    This is implemented as reservoir sampling. The command-line would look
    like:
        
        python poverlap.py fixle wgEncodeRegTfbsClusteredV2.bed.gz CTCF Pol2 --n 10

    For wgEncodeRegTfbsClusteredV2.bed.gz, which contains putative TFBS for
    50+ transcription factors. The output looks like this;

        > observed number of overlaps: 31572
        > shuffle command: bedtools intersect -u -a atype.bed -b <(python poverlap.py bed-sample otypes.bed --n 76623) | wc -l
        > simulated overlap mean: 29895.4
        > simulated p-value: 0
        [30147, 29782, 29835, 30023, 29943, 29844, 29567, 30286, 29868, 29659]

    where, in this case, the shuffle command is repeated 10 times. We see that the
    observed overlap (31752) is higher than any of the simulated. So, the simulated
    p-value is 0.
    
Installation
============

The script is stand-alone, dependencies can be installed with

```Shell
    pip install toolshed commandr
```

Usage
=====

    poverlap is the main function that parallelizes testing overlap between `a`
    and `b`. It performs `n` shufflings and compares the observed number of
    lines in the intersection to the simulated intersections to generate a
    p-value.
    When using shuffle_distance, `exclude`, `include` and `chrom` are ignored.
    Args that are not explicitly part of BEDTools are explained below, e.g. to
    find intervals that are within a given distance, rather than fully
    overlapping, one can set overlap_distance to > 0.
    To shuffle intervals within a certain distance of their current location,
    use shuffle_distance to retain the local structure.

    Arguments:
        a - first bed file
        b - second bed file
        genome - genome file
        n - number of shuffles
        chrom - shuffle within chromosomes
        exclude - optional bed file of regions to exclude
        include - optional bed file of regions to include
        shuffle_both - if set, both a and b are shuffled. normally just b
        overlap_distance - intervals within this distance are overlapping.
        shuffle_distance - shuffle each interval to a random location within
                           this distance of its current location.
