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

# 0      1       2       3       4       5       6       7       8
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  FB      PL      ST      UG      VS

# VCF variant line
class VCFvariant(list):
    def __init__(self,line):
        self.line = line
        try:
            f = line.split()
            self.CHROM = f[0]
            self.POS = f[1]
            self.ID = f[2]
            self.REF = f[3]
            self.ALT = f[4].split(',')
        except:
            print '>>>\n', line, '\n<<<'
            raise
        return

    def __str__(self):
        return "< VCFvariant "+self.CHROM+":"+self.POS+' '+self.REF+"|"+','.join(self.ALT)+" >"

    def locus(self):
        '''returns interval of variant'''
        return ( self.CHROM, int(self.POS), int(self.POS)+len(self.REF) )


class VCF(object):  # reads the whole file!
    def __init__(self,fh):
        self.header = []
        self.samples = None
        self.entries = []
        for line in fh:
            if line.startswith("#") or len(line.rstrip())==0:
                self.header.append(line)
                if line.startswith("#CHROM"):
                    self.samples = line[1:].split()[9:]  # sample header
            else:
                self.entries.append(VCFvariant(line))
        return

    def __iter__(self):
        for e in self.entries:
            yield e
