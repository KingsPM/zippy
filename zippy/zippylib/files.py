#!/usr/bin/env python

__doc__=="""File parsing classes"""
__author__ = "David Brawand"
__license__ = "MIT"
__version__ = "2.3.4"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import sys
import re
import os
from math import ceil
from collections import Counter, defaultdict
from hashlib import sha1
from . import ConfigError
from .interval import *
from urllib import quote, unquote

'''GenePred parser with automatic segment numbering and tiling'''
class GenePred(IntervalList):
    def __init__(self,fh,getgenes=None,interval=None,overlap=None,flank=0,combine=True,noncoding=False):
        IntervalList.__init__(self, [], source='GenePred')
        counter = Counter()
        intervalindex = defaultdict(list)
        # read exons per gene
        genes = defaultdict(list)
        for line in fh:
            if line.startswith("#"):
                continue
            else:
                # create gene and add exons
                f = line.split()
                assert f[3] in ['+','-']
                if getgenes and (f[12] not in getgenes or int(f[6])==int(f[7])) and not noncoding:  # ignore non-coding transcripts
                    continue
                # coding / noncoding
                geneStart = int(f[4]) if noncoding else int(f[6])
                geneEnd = int(f[5]) if noncoding else int(f[7])
                reverse = f[3].startswith('-')
                gene = Interval(f[2],geneStart,geneEnd,f[12],reverse)
                # parse exons
                for e in zip(f[9].split(','),f[10].split(',')):
                    try:
                        map(int,e)
                    except:
                        continue
                    if int(e[1]) < geneStart or geneEnd < int(e[0]):
                        continue  # noncoding
                    try:
                        exonStart = int(e[0]) if noncoding else max(geneStart,int(e[0]))
                        exonEnd = int(e[1]) if noncoding else min(geneEnd,int(e[1]))
                        gene.addSubintervals([Interval(f[2],exonStart,exonEnd,f[12],reverse)])
                    except ValueError:
                        pass
                    except:
                        raise
                # find appropriate gene (same name, and overlapping)
                ovpgenes = [ g for g in genes[gene.name] if gene.overlap(g) ]
                if ovpgenes:
                    try:
                        assert len(ovpgenes) == 1
                    except:
                        # MERGE 2 GENES (there were non-overlapping transcripts in same gene locus!)
                        for i in range(1,len(ovpgenes)):
                            ovpgenes[0].merge(ovpgenes[i],subintervals=True)
                            genes[gene.name].remove(ovpgenes[i])  # remove merged
                    ovpgenes[0].addSubintervals(gene.subintervals)  # add exons from other transcript/gene
                    ovpgenes[0].flattenSubintervals()  # flatten intervals IntervalList
                else:
                    # add new
                    genes[gene.name].append(gene)
        # name metaexons and combine if small enough
        for genename, genelist in genes.items():
            for g in genelist:
                if combine:
                    # iteratively combine closest exons
                    combinedExons = [ [x] for x in sorted(g.subintervals) ]
                    while True:
                        # get distances
                        distances = [ max([ x.chromEnd for x in combinedExons[i]]) - \
                            min([ x.chromStart for x in combinedExons[i-1]]) \
                            for i in range(1,len(combinedExons)) ]
                        # combine smallest distance
                        if any([ d < interval for d in distances ]):
                            smallestIndex = distances.index(min(distances))
                            recombinedExons = []
                            for i,e in enumerate(combinedExons):
                                if i > smallestIndex:  # merge with previous and add remainder
                                    recombinedExons[-1] += e
                                    if i+1 < len(combinedExons):
                                        recombinedExons += combinedExons[i+1:]
                                    break
                                else:
                                    recombinedExons.append(e)
                            combinedExons = recombinedExons
                        else:
                            break
                    # add exon numbers
                    i = 0
                    for e in combinedExons:
                        # get exons
                        ii = range(i,i+len(e))
                        exonNumbers = [ len(g.subintervals) - x for x in ii ] if g.strand < 0 else [ x+1 for x in ii ]
                        if len(e)>1:  # combine exons
                            for j in range(1,len(e)):
                                e[0].merge(e[j])
                        e[0].name += '_{}'.format('+'.join(map(str,sorted(exonNumbers))))
                        intervalindex[e[0].name].append(e[0])
                        i += len(e)
                else:
                    # append exon number
                    for i, e in enumerate(sorted(g.subintervals)):
                        exonNumber = len(g.subintervals) - i if g.strand < 0 else i + 1
                        e.name += '_{}'.format(str(exonNumber))
                        intervalindex[e.name].append(e)
        # split interval if necessary
        for ivs in intervalindex.values():
            for iv in ivs:
                if interval and overlap and interval < len(iv):
                    assert '+' not in iv.name  # paranoia
                    self += iv.tile(interval,overlap,len(f)>3)  # name with suffix if named interval
                else:
                    self += [ iv ]
        # add flanks
        for e in self:
            e.extend(flank)
        return

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
    def __init__(self,fh,flank=0,delim='\t',db=None):
        IntervalList.__init__(self, [], source='VCF')
        self.header = []
        self.samples = []
        self.data = {}
        self.missedgenes = set()
        commentcount = 0
        for i, line in enumerate(fh):
            if line.startswith('#'):
                commentcount += 1
            elif i-commentcount == 0:
                self.header = line.rstrip().split(delim)
                self.data = { h:[] for h in self.header }
            elif re.match('^\s+$',line):
                pass  # tabs/space only line
            else:
                # parse fields
                try:
                    f = line.rstrip().split(delim)
                    row = dict(zip(self.header,f))
                    for k,v in row.items():
                        try:
                            self.data[k].append(v)
                        except:
                            raise Exception('UnknownColumn')
                except:
                    print >> sys.stderr, line
                    print >> sys.stderr, row
                    raise
                # build variant/test description
                if 'primers' in row.keys():  # sample, primer list
                    assert db  # must have database handle to xtract targets
                    pairnames = map(lambda x: x.strip(), row['primers'].split(','))
                    for pp in pairnames:
                        # get pair(s)
                        pairs = db.query(pp)
                        if not pairs:
                            self.missedgenes.add(row['geneID'])
                        # create interval
                        for p in pairs:
                            t = p.sequencingTarget()
                            # (gene,tx,exon,hgvs/pos,zyg)
                            # vd = [ row['geneID'], '', '', '{}:{}-{}'.format(t[0],t[1],t[2]), 'unknown' ]
                            vd = [ row['geneID'], '', '', '', '' ]
                            iv = Interval(t[0],t[1],t[2],name=quote(','.join(vd)),sample=row['sampleID'])
                            iv.extend(-flank)  # shrink search interval
                            self.append(iv)
                else:
                    chrom = row['chromosome'][3:] if row['chromosome'].startswith('chr') else row['chromosome']
                    # parse variant name
                    variantDescription = [ row['geneID'] ]
                    if '-' in row['position']:  # interval (gene,chrom,exon,hgvs/pos,zyg)
                        chromStart, chromEnd = map(int,row['position'].split('-'))
                        variantDescription += [ row['chromosome'] ]
                    else:  # variant (gene,tx,exon,hgvs/pos,zyg)
                        if 'HGVS_c' in row.keys():
                            chromStart, chromEnd = int(row['position']), int(row['position'])+hgvsLength(row['HGVS_c'])
                        elif 'ALT' in row.keys() and 'REF' in row.keys():
                            chromStart, chromEnd = int(row['position']), int(row['position'])+max(map(len,[row['REF'],row['ALT']]))
                        else:
                            raise Exception('UnkownVariantLength')
                        if 'transcriptID' in row.keys():
                            variantDescription += [ row['transcriptID'] ]
                    if 'rank' in row.keys() and '/' in row['rank']:
                        variantDescription += [ 'exon'+row['rank'].split('/')[0] ] # exonnumber
                    variantDescription += [ row['HGVS_c'] if 'HGVS_c' in row.keys() and row['HGVS_c'] else row['position'] ]  # HGVS
                    variantDescription += [ ':'.join([ row[k] for k in sorted(row.keys()) if k.startswith('GT') ]) ]  # zygosity
                    iv = Interval(chrom,chromStart,chromEnd,name=quote(','.join(variantDescription)),sample=row['sampleID'])
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
            if targets.endswith('vcf'):  # VCF files (strand 0)
                intervals = VCF(fh,**tiling)
            elif targets.endswith('bed'):  # BED files (BED3 with automatic names)
                intervals = BED(fh,**tiling)
            elif targets.endswith('txt') or targets.lower().endswith('genepred'):  # GenePred files (uses gene name if unique)
                intervals = GenePred(fh,**tiling)
            else:
                raise Exception('UnknownFileExtension')
    elif re.match('\w+:\d+-\d+',targets):  # single interval, no tiling
        m = re.match('(\w+):(\d+)-(\d+):?([+-])?',targets)
        rev = None if m.group(4) is None else True if m.group(4) == '-' else False
        intervals = [ Interval(m.group(1),m.group(2),m.group(3),reverse=rev) ]
    else:
        raise Exception('FileNotFound')
    return intervals


'''readBatch: read file from SNPpy result output'''
def readBatch(fi,tiling,database=None):
    try:
        assert os.path.isfile(fi)
    except AssertionError:
        print >> sys.stderr, "ERROR: Not a readable file (%s)" % fi
        raise
    with open(fi) as fh:
        intervals = SNPpy(fh,flank=tiling['flank'],db=database)
    sampleVariants = {}
    for iv in intervals:
        try:
            sampleVariants[iv.sample].append(iv)
        except KeyError:
            sampleVariants[iv.sample] = IntervalList([iv],source='SNPpy')
        except:
            raise
    return sampleVariants, \
        sorted(list(set(intervals.data['geneID']))), \
        intervals.missedgenes


'''return length of variant from hgvs.c notation'''
def hgvsLength(hgvs,default=10):
    try:
        m = re.match('c\.-?\d+.+(>|ins|del|dup)(\w+)$',hgvs)
        assert m
    except:
        try:
            l = int(hgvs)
        except:
            print >> sys.stderr, "WARNING: could not find length of variant (%s), assuming %s" % (hgvs,str(default))
            return default
        else:
            return l
    else:
        return len(m.group(2))


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
