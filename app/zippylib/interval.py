#!/usr/bin/env python

__doc__=="""
########################
# interval class/list  #
########################
"""
__author__ = "David Brawand"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import sys
from math import ceil

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
    
