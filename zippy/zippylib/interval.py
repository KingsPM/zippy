#!/usr/bin/env python

__doc__=="""Interval Lists"""
__author__ = "David Brawand"
__license__ = "MIT"
__version__ = "2.3.4"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import sys
from math import ceil

class Interval(object):
    def __init__(self,chrom,chromStart,chromEnd,name=None,reverse=None,sample=None):
        self.chrom = chrom
        self.chromStart = int(chromStart)
        self.chromEnd = int(chromEnd)
        assert self.chromStart <= self.chromEnd  # make sure its on the forward genomic strand
        self.name = name if name else chrom+':'+str(chromStart)+'-'+str(chromEnd)
        self.strand = 0 if reverse is None else -1 if reverse else 1
        self.sample = sample
        self.subintervals = IntervalList([])
        return

    def midpoint(self):
        return int(self.chromStart + (self.chromEnd - self.chromStart)/2.0)

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

    def __repr__(self):
        return "<Interval ("+self.name+") "+self.chrom+":"+str(self.chromStart)+'-'+str(self.chromEnd)+ \
            " ["+str(self.strand)+"] len:"+str(len(self))+">"

    def __str__(self):
        return "\t".join(map(str,[self.chrom, self.chromStart, self.chromEnd, self.name]))

    def tile(self,i,o,suffix=True):  # interval, overlap
        splitintervals = int(ceil( (len(self)-o) / float(i-o) ))  # interval number
        optimalsize = int(ceil( (len(self) + splitintervals*o - o) / float(splitintervals) ))  # optimal interval size
        # get tile spans (and number of exons)
        tilespan = []
        for n,tilestart in enumerate(range(self.chromStart, self.chromEnd, optimalsize-o)):
            tileend = min(tilestart+optimalsize, self.chromEnd)
            tilespan.append((tilestart,tileend))
            if tileend == self.chromEnd:
                break
        tiles = []
        for n,t in enumerate(tilespan):
            tilenumber = len(tilespan)-n if self.strand < 0 else n+1
            tiles.append(Interval(self.chrom,t[0],t[1],self.name+'_'+str(tilenumber) if suffix else None, self.strand < 0))
        return tiles

    def extend(self,flank):
        self.chromStart = self.chromStart-flank if flank <= self.chromStart else 0
        self.chromEnd = self.chromEnd+flank
        return self

    def overlap(self,other):  # also returnd bookended
        return self.chrom == other.chrom and \
            not (other.chromEnd < self.chromStart or other.chromStart > self.chromEnd)

    def merge(self,other,subintervals=False):
        if self.chrom == other.chrom and self.strand == other.strand:
            self.chromStart = other.chromStart if other.chromStart < self.chromStart else self.chromStart
            self.chromEnd = other.chromEnd if other.chromEnd > self.chromEnd else self.chromEnd
            self.name = self.name if other.name == self.name else self.name + '_' + other.name
            if subintervals and (self.subintervals or other.subintervals):
                self.subintervals += other.subintervals
                self.flattenSubintervals()

    def addSubintervals(self,add):
        for e in add:
            if e.chromStart < self.chromStart:
                self.chromStart = e.chromStart
            if e.chromEnd > self.chromEnd:
                self.chromEnd = e.chromEnd
            self.subintervals.append(e)
        self.subintervals.sort()

    def flattenSubintervals(self):
        if self.subintervals:
            self.subintervals.sort()
            merged = [ self.subintervals[0] ]
            for i in range(1,len(self.subintervals)):
                if merged[-1].overlap(self.subintervals[i]):
                    merged[-1].merge(self.subintervals[i])
                else:
                    merged.append(self.subintervals[i])
            self.subintervals = IntervalList(merged)


'''list of intervals'''
class IntervalList(list):
    def __init__(self,elements,source=None):
        list.__init__(self, elements)
        self.source = source  # source of intervals

    def __str__(self):
        return "<IntervalList (%s) %d elements> " % (self.source, len(self))

    def __repr__(self):
        return "<IntervalList (%s) %d elements> " % (self.source, len(self))
