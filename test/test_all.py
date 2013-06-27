from poverlap import poverlap, fixle
import json
import sys

def mymetric(fh):
    val = 0
    for line in fh:
        toks = line.split("\t")
        val += (int(toks[2]) - int(toks[1]))
    return val

def check_attributes(res):
    d = json.loads(res)
    assert len(d) > 0
    d = d.itervalues().next()
    for k in ("sims", "metric", "observed", "shuffle_cmd"):
        assert k in d, (k, d)

def test_python_metric():

    res = poverlap('test/data/a.bed', 'test/data/b.bed', 'data/hg19.genome', metric=mymetric, n=20)
    assert isinstance(res, basestring)
    yield check_attributes, res

def test_string_metric():
    res = poverlap('test/data/a.bed', 'test/data/b.bed', 'data/hg19.genome',
            metric='wc -l', n=20)
    assert isinstance(res, basestring)
    yield check_attributes, res

def test_map_fn():
    res = poverlap('test/data/a.bed', 'test/data/b.bed', 'data/hg19.genome',
            metric='wc -l', n=20, ncpus=map)
    assert isinstance(res, basestring)
    yield check_attributes, res
    from itertools import imap
    res = poverlap('test/data/a.bed', 'test/data/b.bed', 'data/hg19.genome',
            metric='wc -l', n=20, ncpus=imap)
    assert isinstance(res, basestring)
    yield check_attributes, res

def test_cpu_count():

    res = poverlap('test/data/a.bed', 'test/data/b.bed', 'data/hg19.genome',
            metric='wc -l', n=20, ncpus=2)
    assert isinstance(res, basestring)
    yield check_attributes, res

def test_multi_metric():
    res = poverlap('test/data/a.bed', 'test/data/b.bed', 'data/hg19.genome',
            metric=('wc -l', mymetric), n=20, ncpus=2)
    assert isinstance(res, basestring)
    d = json.loads(res)
    assert len(d) == 2, d




