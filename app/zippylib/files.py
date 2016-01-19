#!/usr/bin/env python

__doc__=="""
######################
# simple VCF parser  #
######################
"""
__author__ = "David Brawand"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import sys
import re
import os
from math import ceil
from collections import Counter, defaultdict
from hashlib import sha1
from zippylib import ConfigError

class Interval(object):
    def __init__(self,chrom,chromStart,chromEnd,name=None,reverse=False,sample=None):
        self.chrom = chrom
        self.chromStart = int(chromStart)
        self.chromEnd = int(chromEnd)
        self.name = name if name else chrom+':'+str(chromStart)+'-'+str(chromEnd)
        self.strand = -1 if reverse else 1
        self.sample = sample
        return

    def locus(self):
        '''returns interval of variant'''
        return ( self.chrom, self.chromStart, self.chromEnd )

    def __hash__(self):
        return hash(str(self))

    def __len__(self):
        return self.chromEnd - self.chromStart

    def __eq__(self,other):
        return hash(self) == hash(other)

    def __lt__(self,other):
        return (self.chrom, self.chromStart, self.chromEnd) < (other.chrom, other.chromStart, other.chromEnd)

    def __str__(self):
        return "<Interval ("+self.name+") "+self.chrom+":"+str(self.chromStart)+'-'+str(self.chromEnd)+ \
            " ["+str(self.strand)+"] len:"+str(len(self))+">"

    def tile(self,i,o,suffix=True):  # interval, overlap
        splitintervals = int(ceil( (len(self)-o) / float(i-o) ))  # interval number
        optimalsize = int(ceil( (len(self) + splitintervals*o - o) / float(splitintervals) ))  # optimal interval size
        tiles = []
        for n,tilestart in enumerate(range(self.chromStart, self.chromEnd, optimalsize-o)):
            tileend = min(tilestart+optimalsize, self.chromEnd)
            tiles.append(Interval(self.chrom,tilestart,tileend, self.name+'_'+str(n+1) if suffix else None, self.strand < 0))
            if tileend == self.chromEnd:
                break
        return tiles

    def extend(self,flank):
        self.chromStart = self.chromStart-flank if flank <= self.chromStart else 0
        self.chromEnd = self.chromEnd+flank
        return self

'''list of intervals'''
class IntervalList(list):
    def __init__(self,elements,source=None):
        list.__init__(self, elements)
        self.source = source  # source of intervals

    def __str__(self):
        return "<IntervalList (%s) %d elements> " % (self.source, len(self))


'''bed parser with automatic segment numbering and tiling'''
class BED(IntervalList):
    def __init__(self,fh,interval=None,overlap=None,flank=0):
        IntervalList.__init__(self, [], source='BED')
        counter = Counter()
        intervalindex = defaultdict(list)
        for line in fh:
            if line.startswith("#"):
                continue
            else:
                # create interval
                f = line.split()
                try:
                    if len(f) > 5:  # name/strand
                        iv = Interval(f[0],int(f[1]),int(f[2]),f[3],f[5].startswith('-'))
                    elif len(f) > 3:  # name
                        iv = Interval(f[0],int(f[1]),int(f[2]),f[3])
                    else:  # automatic naming
                        iv = Interval(f[0],int(f[1]),int(f[2]))
                except:
                    print >> sys.stderr, f
                    raise
                intervalindex[iv.name].append(iv)
        # suffix interval names if necessary
        for ivs in intervalindex.values():
            if len(ivs)>1:
                for i,iv in enumerate(ivs):
                    iv.name += '-{:02d}'.format(i+1)
        # split interval if necessary
        for ivs in intervalindex.values():
            for iv in ivs:
                if interval and overlap and interval < len(iv):
                    self += iv.tile(interval,overlap,len(f)>3)  # name with suffix if named interval
                else:
                    self += [ iv ]
        # add flanks
        for e in self:
            e.extend(flank)
        return

'''vcf parser with segment hashing and tiling'''
class VCF(IntervalList):  # no interval tiling as a variant has to be sequenced in onse single run
    def __init__(self,fh,interval=None,overlap=None,flank=0):
        IntervalList.__init__(self, [], source='VCF')
        self.header = []
        self.samples = None
        for line in fh:
            if line.startswith("#") or len(line.rstrip())==0:
                self.header.append(line)
                if line.startswith("#CHROM"):
                    self.samples = line[1:].split()[9:]  # sample header
            else:
                f = line.split()
                iv = Interval(f[0],int(f[1]),int(f[1])+max(map(len,[f[3]]+f[4].split(','))),name=f[2] if f[2]!='.' else None)
                self.append(iv)
        # add flanks and name
        for e in self:
            e.extend(flank)
        return

'''SNPpy result reader'''
class SNPpy(IntervalList):
    def __init__(self,fh,flank=0,delim='\t'):
        IntervalList.__init__(self, [], source='VCF')
        self.header = []
        self.samples = []
        for i, line in enumerate(fh):
            if i == 0:
                self.header = line.split(delim)
            else:
                f = line.split()
                row = dict(zip(self.header,f))
                chrom = row['chromosome'][3:] if row['chromosome'].startswith('chr') else row['chromosome']
                chromStart = int(row['position'])
                chromEnd = chromStart+hgvsLength(row['HGVS_c'])
                variantName = ':'.join([row['geneID'],row['transcriptID'],row['HGVS_c'],row['GT/CONSENSUS']]).replace('>','to')
                iv = Interval(chrom,chromStart,chromEnd,name=variantName,sample=row['sampleID'])
                self.append(iv)
        # add flanks and name
        for e in self:
            e.extend(flank)
        return


'''generic data class with formatted output'''
class Data(object):
    def __init__(self,data,header):
        self.header = header
        self.data = data

    def writefile(self, fi):
        if fi.endswith('.interval') or fi.endswith('.bed'):
            try:
                assert set(['chrom','chromStart','chromEnd','name']).issubset(set(self.header))
                fh = open(fi,'w')
            except AssertionError:
                raise ConfigError('Wrong data format for interval/bed')
            except:
                raise
            else:
                for row in sorted(self.data):
                    d = dict(zip(self.header,row))
                    if fi.endswith('bed'):
                        fh.write('\t'.join(map(str,[d['chrom'], d['chromStart'], d['chromEnd'], d['name']]))+'\n')
                    else:
                        fh.write('\t'.join(map(str,[d['chrom'], d['chromStart'], d['chromEnd'], '+', d['name']]))+'\n')
                fh.close()
        else: # write as is
            fh = sys.stdout if fi == '-' else open(fi,'w')
            print >> fi, self.header
            for d in self.data:
                print >> fi, '\t'.join(map(str,[ d[f] for f in self.header]))
            if fi != '-':
                fh.close()


''' read target intervals from VCF, BED or directly'''
def readTargets(targets,tiling):
    if os.path.isfile(targets):
        with open(targets) as fh:
            if targets.endswith('vcf'):  # VCF files
                intervals = VCF(fh,**tiling)
            elif targets.endswith('bed'):  # BED files (BED3 with automatic names)
                intervals = BED(fh,**tiling)
            else:
                raise Exception('UnknownFileExtension')
    elif re.match('\w+:\d+-\d+',targets):  # single interval
        m = re.match('(\w+):(\d+)-(\d+)',targets)
        intervals = [ Interval(*m.groups()) ]
    else:
        Exception('FileNotFound')
    return intervals


'''readBatch: read file from SNPpy result output'''
def readBatch(fi,tiling):
    try:
        assert os.path.isfile(fi)
    except AssertionError:
        print >> sys.stderr, "ERROR: Not a readable file (%s)" % fi
        raise
    with open(fi) as fh:
        intervals = SNPpy(fh,flank=tiling['flank'])
    sampleVariants = {}
    for iv in intervals:
        try:
            sampleVariants[iv.sample].append(iv)
        except KeyError:
            sampleVariants[iv.sample] = IntervalList([iv],source='SNPpy')
        except:
            raise
    return sampleVariants


'''return length of variant from hgvs.c notation'''
def hgvsLength(hgvs,default=10):
    try:
        assert hgvs.startswith('c.')
    except:
        raise
    try:
        m = re.match('c.\d+(\D+)(>|ins|del|dup)(\w+)$',hgvs)
        assert m
    except:
        print >> sys.stderr, "WARNING: could not find length of variant, assuming %s" % str(default)
        return default
    else:
        return len(m.group(3))


if __name__=="__main__":
    cfg = {
        'interval': 200,
        'overlap': 50,
        'flank': 35
    }
    with open(sys.argv[1]) as fh:
        bed = BED(fh,**cfg)
        for i in bed:
            print i
