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

    a) allows shuffing intervals in `a` to a new location inside the interval
       that contains them in `b`.

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
        
    ./poverlap.py fixle wgEncodeRegTfbsClusteredV2.bed.gz CTCF Pol2 --n 10

For wgEncodeRegTfbsClusteredV2.bed.gz, which contains putative TFBS for
50+ transcription factors. The output looks like this;

    > observed number of overlaps: 31572
    > shuffle command: bedtools intersect -u -a atype.bed -b <(./poverlap.py bed-sample otypes.bed --n 76623) | wc -l
    > simulated overlap mean: 29895.4
    > simulated p-value: 0
    [30147, 29782, 29835, 30023, 29943, 29844, 29567, 30286, 29868, 29659]

where, in this case, the shuffle command is repeated 10 times. We see that the
observed overlap (31752) is higher than any of the simulated. So, the simulated
p-value is 0.

Examples
========
    
shuffle b to within X bases of original location

    ./poverlap.py poverlap --a a.bed --b b.bed \
            --shuffle_loc 10000 \
            > res.shuffle_loc-10kb.txt

shuffle both a and b to within X bases of original location

    ./poverlap.py poverlap --a a.bed --b b.bed \
            --shuffle_both --shuffle_loc 10000 \
            > res.both.shuffle_loc-10kb.txt

we can also shuffle intervals in b to a random location within the interval
that contains it in 'other.bed'. For example, b.bed may be transcription start
sites and we wish to see what things look like if we randomize the location
of the TSS to anywhere within the gene-body.

    ./poverlap.py poverlap --a a.bed --b b.bed \
            --shuffle_loc gene-bodies.bed \
            > res.gene-bodies.txt

It is possible to use Haiminen's method to shuffle to known sites.
Here, a and b are a pair of TFs out of the 50+ in the
wgEncodeRegTfbsClusteredV2.bed.gz

    ./poverlap.py fixle \
        wgEncodeRegTfbsClusteredV2.bed.gz a.bed b.bed > res.fixle.haimenin.txt

simple bedtools shuffle

    ./poverlap.py poverlap --a a.bed --b b.bed \
        -g hg19.genome > res.shuffle-bt.txt

bedtools exclude:

    ./poverlap.py poverlap --a a.bed --b b.bed \
        -g hg19.genome --exclude repressed.encode.bed > res.shuffle_bt.exclude-repressed.txt

bedtools include:

    ./poverlap.py poverlap --a a.bed --b b.bed \
        -g hg19.genome --include hg19.promoters.bed > res.shuffle_bt.include-promoters.txt

bedtools include and exclude:

    ./poverlap.py poverlap --a a.bed --b b.bed \
        -g hg19.genome --include hg19.promoters.bed --exclude repressed.encode.bed \
        > res.shuffle_bt.include-promoters.exclude-repressed.txt

Call intervals within 40kb of another as "overlapping". Note, this is different
from shuffle\_distance which determines where the interval can go.

    ./poverlap.py poverlap --a a.bed --b b.bed \
        --overlap_distance 40000 -g hg19.genome > res.shuffle-bt-within-40kb.txt

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

The best way to understand is to use the script. E.g start with:

    ./poverlap.py poverlap

When using shuffle_loc, `exclude`, `include` and `chrom` are ignored.
Args that are not explicitly part of BEDTools are explained below, e.g. to
find intervals that are within a given distance, rather than fully
overlapping, one can set overlap_distance to > 0.
To shuffle intervals within a certain distance of their current location,
use shuffle_loc to retain the local structure.

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
    shuffle_loc - randomize the start of each interval to a random
                  location within this distance of its current location.
                   or if a BED file is given, then randomize the intervals in
                   `b` to an random location within the containing interval in
                   shuffle_loc
