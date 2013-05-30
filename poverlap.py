
from toolshed import nopen
import sys
from multiprocessing import Pool
import argparse
from commandr import command, Run

def run(cmd):
    return list(nopen("|%s" % cmd.lstrip("|")))[0]

def extend_bed(fin, fout, bases):
    bases /= 2
    with open(fout, 'w') as fh:
        for toks in (l.rstrip("\r\n").split("\t") for l in nopen(fin)):
            toks[1] = max(0, int(toks[1]) - bases)
            toks[2] = max(0, int(toks[2]) + bases)
            if toks[1] > toks[2]: # negative distances
                toks[1] = toks[2] = (toks[1] + toks[2]) / 2
            assert toks[1] <= toks[2]
            print >>fh, "\t".join(map(str, toks))
    return fh.name

@command('fixle') # from Haiminen et al in BMC Bioinformatics 2008, 9:336
def fixle(bed, atype, btype, type_col=4, n=10):
    """
    `bed` may contain, e.g. 20 TFBS as defined by the type in `type_col`
    we keep the rows labeld as `atype` in the same locations, but we randomly
    assign `btype` to any of the remaining rows.
    Arguments:
        bed - BED file with a column that delineates types
        atype - the query type, e.g. Pol2
        btype - the type to be shuffled, e.g. CTCF
        type_col - the column in `bed` the lists the types
        n - number of shuffles
    """
    type_col -= 1
    n_btypes = 0
    pool = Pool(10)
    with open('atype.bed', 'w') as afh, \
            open('otypes.bed', 'w') as ofh, \
            open('btypes.bed', 'w') as bfh:
        for toks in (l.rstrip("\r\n").split("\t") for l in nopen(bed)):
            if toks[type_col] == atype:
                print >> afh, "\t".join(toks)
            else:
                print >> ofh, "\t".join(toks)
                if toks[type_col] == btype:
                    print >>bfh, "\t".join(toks)
                    n_btypes += 1
    assert n_btypes > 0, ("no intervals found for", btype)

    a, b, other = afh.name, bfh.name, ofh.name
    orig_cmd = "bedtools intersect -u -a {a} -b {b} | wc -l".format(**locals())
    observed = int(run(orig_cmd))
    print "> observed number of overlaps: %i" % observed
    script = __file__
    bsample = '<(python {script} bed-sample {other} --n {n_btypes})'.format(**locals())
    shuf_cmd = "bedtools intersect -u -a {a} -b {bsample} | wc -l".format(**locals())
    print "> shuffle command: %s" % shuf_cmd
    sims = [int(x) for x in pool.imap(run, [shuf_cmd] * n)]
    print "> simulated overlap mean: %.1f" % (sum(sims) / float(len(sims)))
    print "> simulated p-value: %.3g" \
            % (sum((s >= observed) for s in sims) / float(len(sims)))

    print sims


@command('bed-sample')
def bed_sample(bed, n=100):
    """
    choose n random lines from a bed file.
    uses reservoir sampling
    Arguments:
        bed - a bed file
        n - number of lines to sample [default: %default]
    """
    n = int(n)
    import random
    lines = []
    for i, line in enumerate(nopen(bed)):
        if i < n:
            lines.append(line)
        else:
            replace_idx = random.randint(0, i)
            if replace_idx < n:
                lines[replace_idx] = line
    print "".join(lines),

@command('distance-shuffle')
def distance_shuffle(bed, dist=500000):
    """
    randomize the location of each interval in `bed` by moving it's
    start location to within `dist` bp of its current location.
    Arguments:
        bed - input bed file
        dist - shuffle intervals to within this distance (+ or -)
    """
    from random import randint
    dist = abs(int(dist))
    for toks in (l.rstrip('\r\n').split('\t') for l in nopen(bed)):
        d = randint(-dist, dist)
        toks[1:3] = [str(max(0, int(loc) + d)) for loc in toks[1:3]]
        print "\t".join(toks)

@command('poverlap')
def poverlap(a, b, genome=None, n=50, chrom=False, exclude=None, include=None,
        shuffle_both=False, overlap_distance=0, shuffle_distance=None):
    """\
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
    """
    pool = Pool(8)

    n = int(n)
    exclude = "" if exclude is None else ("-excl %s" % exclude)
    include = "" if include is None else ("-incl %s" % include)
    chrom = "" if chrom is False else "-chrom"
    if genome is None: assert shuffle_distance

    if exclude: # remove input regions that fall in exclude
        run("bedtools intersect -v -a {a} -b {e} > aa.bed; echo 1"\
                .format(a=a, e=exclude.split()[1]))
        a = "aa.bed"
        run("bedtools intersect -v -a {b} -b {e} > bb.bed; echo 1"\
                .format(b=b, e=exclude.split()[1]))
        b = "bb.bed"

    if include: # remove input regions that DONT fall in exclude
        run("bedtools intersect -u -a {a} -b {i} > aaa.bed; echo 1"\
                .format(a=a, i=include.split()[1]))
        a = "aaa.bed"
        run("bedtools intersect -u -a {b} -b {i} > bbb.bed; echo 1"\
                .format(b=b, i=include.split()[1]))
        b = "bbb.bed"

    if overlap_distance != 0:
        a = extend_bed(a, 'aaaa.bed', overlap_distance)
        b = extend_bed(b, 'bbbb.bed', overlap_distance)

    orig_cmd = "bedtools intersect -u -a {a} -b {b} | wc -l".format(**locals())

    if shuffle_distance is None:
        # use bedtools shuffle
        if shuffle_both:
            a = "<(bedtools shuffle {exclude} {include} -i {a} -g {genome} {chrom})".format(**locals())
        shuf_cmd = ("bedtools intersect -u -a {a} "
                "-b <(bedtools shuffle {exclude} {include} -i {b} -g {genome} {chrom}) | wc -l ".format(**locals()))
    else:
        # use python shuffle ignores --chrom and --genome
        shuffle_distance = int(shuffle_distance)
        script = __file__
        if shuffle_both:
            a = "<(python {script} distance-shuffle {a} --dist {shuffle_distance})".format(**locals())
        shuf_cmd = ("bedtools intersect -u -a {a} "
            "-b <(python {script} distance-shuffle {b} --dist {shuffle_distance})"
            " | wc -l").format(**locals())

    #print "original command: %s" % orig_cmd
    print "> shuffle command: %s" % shuf_cmd
    observed = int(run(orig_cmd))
    print "> observed number of overlaps: %i" % observed
    sims = [int(x) for x in pool.imap(run, [shuf_cmd] * n)]
    print "> simulated overlap mean: %.1f" % (sum(sims) / float(len(sims)))
    print "> simulated p-value: %.3g" % (sum((s >= observed) for s in sims) / float(len(sims)))
    print ">", sims

if __name__ == "__main__":
    res = Run()
