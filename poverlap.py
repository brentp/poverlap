#!/usr/bin/env python

import os
import sys
import json
from commandr import command, Run
from toolshed import nopen, reader
from multiprocessing import cpu_count
from multiprocessing.dummy import Pool
from tempfile import mktemp as _mktemp

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

NCPUS = cpu_count()
if NCPUS > 4: NCPUS -= 1

SEP = "&*#Z"


def mktemp(*args, **kwargs):
    def rm(f):
        try: os.unlink(f)
        except OSError: pass

    if not 'suffix' in kwargs: kwargs['suffix'] = ".bed"
    f = _mktemp(*args, **kwargs)
    atexit.register(rm, f)
    return f


def run(cmd):
    proc = nopen("|%s" % cmd.lstrip("|"), mode=None)
    ret = proc.stdout.next()
    check_proc(proc, cmd)
    return ret


def check_proc(proc, cmd=''):
    err = proc.stderr.read().strip()
    proc.terminate()
    if proc.returncode not in (0, None):
        sys.stderr.write("%s\n%s\n\%s" % (cmd, err, proc.returncode))
        raise Exception(err)
    if err: sys.stderr.write(err)
    proc.stdout.close(); proc.stderr.close()
    del proc


def run_metric(cmd, metric=None):
    """
    Metric can be a string, e.g. "wc -l" or a python callable that consumes
    lines of input and returns a single value.
    e.g.

    def mymetric(fh):
        val = 0
        for line in fh:
            val += float(line.split("\t")[4])
        return val

    The lines sent to the metric function will be the result of bedtools
    intersect -wa -- so that both the -a and -b intervals will be present
    in each line.
    """

    if metric is None:
        cmd, metric = cmd
    if isinstance(metric, basestring):
        return float(run("%s | %s" % (cmd, metric)))
    else:
        proc = nopen("|%s" % cmd, mode=None)
        res = metric(proc.stdout)
        check_proc(proc, cmd)
        assert isinstance(res, (int, float))
        return res


def extend_bed(fin, fout, bases):
    # `bedtools slop`

    # we're extending both a.bed and b.bed by this distance
    # so divide by 2.
    bases /= 2
    with nopen(fout, 'w') as fh:
        for toks in (l.rstrip("\r\n").split("\t") for l in nopen(fin)):
            toks[1] = max(0, int(toks[1]) - bases)
            toks[2] = max(0, int(toks[2]) + bases)
            if toks[1] > toks[2]:  # negative distances
                toks[1] = toks[2] = (toks[1] + toks[2]) / 2
            assert toks[1] <= toks[2]
            print >>fh, "\t".join(map(str, toks))
    return fh.name


@command('fixle')
def fixle(bed, atype, btype, type_col=4, metric='wc -l', n=100, ncpus=-1):
    """\
    From Haiminen et al in BMC Bioinformatics 2008, 9:336 (and R's `cooccur`).
    `bed` may contain, e.g. 20 TFBS as defined by the type in `type_col`.
    We keep the rows labeled as `atype` in the same locations, but we randomly
    assign `btype` to any of the remaining rows.

    Arguments:
        bed - BED file with a column that delineates types
        atype - the query type, e.g. Pol2
        btype - the type to be shuffled, e.g. CTCF
        type_col - the column in `bed` the lists the types
        n - number of shuffles
        metric - a string that indicates a program that consumes BED intervals
        ncpus - number cpus to use -- if a callable does the parallelization
                use, e.g. Pool(5).map or Ipython Client[:].map
    """
    type_col -= 1
    n_btypes = 0
    pmap = get_pmap(ncpus)
    with nopen(mktemp(), 'w') as afh, \
            nopen(mktemp(), 'w') as ofh, \
            nopen(mktemp(), 'w') as bfh:
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
    orig_cmd = "bedtools intersect -wa -a {a} -b {b}".format(**locals())
    script = __file__
    bsample = '<(python {script} bed-sample {other} --n {n_btypes})'.format(**locals())
    shuf_cmd = "bedtools intersect -wa -a {a} -b {bsample}".format(**locals())
    return json.dumps(gen_results(orig_cmd, metric, pmap, n, shuf_cmd))


@command('bed-sample')
def bed_sample(bed, n=100):
    """\
    Choose n random lines from a bed file. Uses reservoir sampling.

    Arguments:
        bed - a bed file
        n - number of lines to sample
    """
    n, lines = int(n), []
    from random import randint
    with nopen(bed) as fh:
        for i, line in enumerate(nopen(fh)):
            if i < n:
                lines.append(line)
            else:
                replace_idx = randint(0, i)
                if replace_idx < n:
                    lines[replace_idx] = line
        print "".join(lines),


@command('local-shuffle')
def local_shuffle(bed, loc='500000'):
    """
    Randomize the location of each interval in `bed` by moving its
    start location to within `loc` bp of its current location or to
    its containing interval in `loc`.

    Arguments:
        bed - input bed file
        loc - shuffle intervals to within this distance (+ or -).
               If not an integer, then this should be a BED file containing
               regions such that each interval in `bed` is shuffled within
               its containing interval in `loc`
    """
    from random import randint
    if str(loc).isdigit():
        dist = abs(int(loc))
        with nopen(bed) as fh:
            for toks in (l.rstrip('\r\n').split('\t') for l in fh):
                d = randint(-dist, dist)
                toks[1:3] = [str(max(0, int(bloc) + d)) for bloc in toks[1:3]]
                print "\t".join(toks)
    else:
        # we are using dist as the windows within which to shuffle
        assert os.path.exists(loc)
        bed4 = mktemp()
        with open(bed4, 'w') as fh:
            # this step is so we don't have to track the number of columns in A
            for toks in reader(bed, header=False):
                fh.write("%s\t%s\n" % ("\t".join(toks[:3]), SEP.join(toks)))

        missing = 0
        # we first find the b-interval that contains each a-interval by
        # using bedtools intersect
        for toks in reader("|bedtools intersect -wao -a {bed4} -b {loc}"
                           .format(**locals()), header=False):
            ajoin = toks[:4]
            a = ajoin[3].split(SEP)  # extract the full interval
            b = toks[4:]

            if int(b[-1]) == 0:
                missing += 1
                continue
            assert a[0] == b[0], ('chroms dont match', a, b)

            alen = int(a[2]) - int(a[1])
            # doesn't care if the new interval is completely contained in b
            astart = randint(int(b[1]), int(b[2]))

            # subtract half the time.
            aend = (astart - alen) if randint(0, 1) == 0 and astart > alen \
                else (astart + alen)

            a[1], a[2] = map(str, (astart, aend) if astart < aend
                             else (aend, astart))

            print "\t".join(a)
        if missing > 0:
            print >> sys.stderr, ("found {missing} intervals in {bed} that "
                                  " were not contained in {loc}"
                                  .format(**locals()))


def zclude(bed, other, exclude=True):
    """
    Include or exclude intervals from bed that overlap other.
    If exclude is True:
        new = bedtools intersect -v -a bed -o other
    """
    if other is None: return bed
    n_orig = sum(1 for _ in nopen(bed))
    tmp = mktemp()
    if exclude:
        run("bedtools intersect -v -a {bed} -b {other} > {tmp}; echo 1"
            .format(**locals()))
    else:
        run("bedtools intersect -u -a {bed} -b {other} > {tmp}; echo 1"
            .format(**locals()))
    n_after = sum(1 for _ in nopen(tmp))
    clude = "exclud" if exclude else "includ"
    pct = 100 * float(n_orig - n_after) / n_orig
    print >>sys.stderr, ("reduced {bed} from {n_orig} to {n_after} "
             "{pct:.3f}% by {clude}ing {other}").format(**locals())
    return tmp


def get_pmap(ncpus):
    if ncpus in ('1', 1, None):
        pmap = map
    elif isinstance(ncpus, (basestring, int)):
        ncpus = int(ncpus)
        if ncpus == -1: ncpus = cpu_count()
        pool = Pool(ncpus)

        ############################################################
        # this block seems to be necessary to avoid errors at exit #
        ############################################################
        import atexit
        def term():
            try: pool.terminate()
            except: pass
        atexit.register(term)
        ############################################################

        pmap = pool.imap
    else:
        pmap = ncpus
        assert hasattr(pmap, "__call__"), pmap
    return pmap


@command('poverlap')
def poverlap(a, b, genome=None, metric='wc -l', n=100, chrom=False,
             exclude=None, include=None, shuffle_both=False,
             overlap_distance=0, shuffle_loc=None, ncpus=-1):
    """\
    poverlap is the main function that parallelizes testing overlap between `a`
    and `b`. It performs `n` shufflings and compares the observed number of
    lines in the intersection to the simulated intersections to generate a
    p-value.

    When using shuffle_loc, `exclude`, `include` and `chrom` are ignored.
    Args that are not explicitly part of BEDTools are explained below, e.g. to
    find intervals that are within a given distance, rather than fully
    overlapping, one can set overlap_distance to > 0. To shuffle intervals
    within a certain distance of their current location, or to keep them
    inside a set of intervals, use shuffle_loc to retain the local structure.

    Arguments:
        a - first bed file
        b - second bed file
        genome - genome file
        metric - a string that indicates a program that consumes BED intervals
                 from STDIN and outputs a single, numerical value upon
                 completion. Default is 'wc -l'
        n - number of shuffles
        chrom - shuffle within chromosomes
        exclude - optional bed file of regions to exclude
        include - optional bed file of regions to include
        shuffle_both - if set, both A and B are shuffled. Default is B only.
        overlap_distance - intervals within this distance are overlapping
        shuffle_loc - shuffle each interval to a random location within this
                      distance of its current location. If not an integer,
                      then this should be a BED file containing regions such
                      that each interval in `bed` is shuffled within its
                      containing interval in `shuffle_loc`.
        ncpus - number cpus to use -- if a callable does the parallelization
                use, e.g. Pool(5).map or Ipython Client[:].map
    """
    pmap = get_pmap(ncpus)

    assert os.path.exists(genome), (genome, "not available")

    n = int(n)
    chrom = "" if chrom is False else "-chrom"
    if genome is None: assert shuffle_loc

    # limit exclude and then to include
    a = zclude(zclude(a, exclude, True), include, False)
    b = zclude(zclude(b, exclude, True), include, False)

    exclude = "" if exclude is None else ("-excl %s" % exclude)
    include = "" if include is None else ("-incl %s" % include)

    if overlap_distance != 0:
        a = extend_bed(a, mktemp(), overlap_distance)
        b = extend_bed(b, mktemp(), overlap_distance)

    orig_cmd = "bedtools intersect -wa -a {a} -b {b}".format(**locals())

    if shuffle_loc is None:
        # use bedtools shuffle
        if shuffle_both:
            a = ("<(bedtools shuffle {exclude} {include} -i {a} -g {genome} "
                 "{chrom})".format(**locals()))
        shuf_cmd = ("bedtools intersect -wa -a {a} -b "
                    "<(bedtools shuffle {exclude} {include} -i {b} -g {genome}"
                    " {chrom})".format(**locals()))
    else:
        # use python shuffle ignores --chrom and --genome
        script = __file__
        if shuffle_both:
            a = "<(python {script} local-shuffle {a} --loc {shuffle_loc})"\
                .format(**locals())
        shuf_cmd = ("bedtools intersect -wa -a {a} -b "
                    "<(python {script} local-shuffle {b} --loc {shuffle_loc})"
                    ).format(**locals())

    return json.dumps(gen_results(orig_cmd, metric, pmap, n, shuf_cmd))


def gen_results(orig_cmd, metric, pmap, n, shuf_cmd=None):
    if not isinstance(metric, (tuple, list)):
        metric = [metric]
    full_res = {}
    for met in metric:
        observed = run_metric(orig_cmd, met)
        res = {"observed": observed, "shuffle_cmd": shuf_cmd}
        sims = [int(x) for x in pmap(run_metric, [(shuf_cmd, met)] * n)]
        res['metric'] = repr(met)
        res['simulated mean metric'] = (sum(sims) / float(len(sims)))
        res['simulated_p'] = \
            (sum((s >= observed) for s in sims) / float(len(sims)))
        res['sims'] = sims
        full_res[repr(met)] = res
    return full_res


def main():
    res = Run()


if __name__ == "__main__":
    main()
