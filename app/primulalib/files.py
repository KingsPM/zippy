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
from math import ceil
from collections import Counter
from hashlib import md5

class Interval(object):
    def __init__(self,chrom,chromStart,chromEnd,name='',reverse=False):
        self.chrom = chrom
        self.chromStart = chromStart
        self.chromEnd = chromEnd
        self.name = name
        self.strand = -1 if reverse else 1
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

    def __str__(self):
        return "<Interval ("+self.name+") "+self.chrom+":"+str(self.chromStart)+'-'+str(self.chromEnd)+ \
            " ["+str(self.strand)+"] len:"+str(len(self))+">"

    def tile(self,i,o):
        splitintervals = int(ceil((len(self)+o) / float(i)))  #integer division
        optimalsize = int(ceil( (len(self) + splitintervals*o) / float(1 + splitintervals)))
        tiles = []
        for tilestart in range(self.chromStart, self.chromEnd, optimalsize-o):
            tileend = min(tilestart+optimalsize, self.chromEnd)
            tiles.append(Interval(self.chrom,tilestart,tileend, self.name, self.strand < 0))
            if tileend == self.chromEnd:
                break
        return tiles

    def extend(self,flank):
        self.chromStart = self.chromStart-flank if flank <= self.chromStart else 0
        self.chromEnd = self.chromEnd+flank
        return self

'''bed parser with automatic segment numbering and tiling'''
class BED(object):
    def __init__(self,fh,interval=None,overlap=None,flank=0):
        self.entries = []
        counter = Counter()
        tiles = []
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
                    else:  # anonymous interval
                        iv = Interval(f[0],int(f[1]),int(f[2]))
                except:
                    print >> sys.stderr, f
                    raise
                # tile interval
                if interval and overlap and interval < len(iv):
                    tiles += iv.tile(interval,overlap)
                else:
                    tiles += [ iv ]
        # count opposite strand tiles
        for t in tiles:
            if t.strand < 0:
                counter[t.name] += 1
            else:
                counter[t.name] = 1
        # append number to interval names
        for t in tiles:
            name = t.name
            t.name += '_'+str(counter[name])
            counter[name] += t.strand
        # add flanks
        for tile in tiles:
            tile.extend(flank)
        # store
        self.entries = tiles
        return

    def __iter__(self):
        for e in self.entries:
            yield e

'''vcf parser with segment hashing and tiling'''
class VCF(object):  # reads the whole file!
    def __init__(self,fh,interval=None,overlap=None,flank=0):
        self.header = []
        self.samples = None
        self.entries = []
        for line in fh:
            if line.startswith("#") or len(line.rstrip())==0:
                self.header.append(line)
                if line.startswith("#CHROM"):
                    self.samples = line[1:].split()[9:]  # sample header
            else:
                f = line.split()
                iv = Interval(f[0],int(f[1]),int(f[1])+len(f[3]))
                if interval and overlap and interval < len(iv):
                    for tile in iv.tile(interval,overlap):
                        self.entries.append(tile)
                else:
                    self.entries.append(iv)
        # add flanks and name
        for e in self.entries:
            c.extend(flank)
            e.name = 'X'+md5((e.chrom,e.chromStart,e.chromEnd)).hexdigest()[:8]
        return

    def __iter__(self):
        for e in self.entries:
            yield e


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
