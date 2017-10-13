#!/usr/bin/env python
# -*- coding: utf-8 -*-

__doc__=="""Primer3 Classes"""
__author__ = "David Brawand"
__license__ = "MIT"
__version__ = "2.3.4"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import sys, os, re, datetime
from hashlib import md5, sha1
import primer3
import pysam
import subprocess
from collections import defaultdict, OrderedDict, Counter
from .interval import Interval
from string import maketrans
from urllib import unquote
revcmp = maketrans('ACGTNacgtn','TGCANtgcan')

'''returns common prefix (substring)'''
def commonPrefix(left,right,stripchars='-_ ',commonlength=3):
    if left and right:
        matchingPositions = [ i+1 for i,j in enumerate([ i for i, x in enumerate(zip(left,right)) if len(set(x)) == 1]) if i==j]
        if matchingPositions and max(matchingPositions) >= commonlength:
            return left[:max(matchingPositions)].rstrip(stripchars)
    return None

'''return -1,0,1'''
def parsePrimerName(x):
    fwd_suffix, rev_suffix = ['f','fwd','5\'','left'], ['r','rev','3\'','right']
    if x[0:2] in ['3\'','5\'']:  # prefix case
        if x.startswith('5'):
            return (x[2:], 1)
        elif x.startswith('3'):
            return (x[2:], -1)
    else:  # suffix case
        pre, suf = x[:x.replace('-','_').rfind('_')], x[x.replace('-','_').rfind('_')+1:].lower()
        if suf in fwd_suffix or re.match(r'f\d+',suf):
            return (pre, 1)
        elif suf in rev_suffix or re.match(r'r\d+',suf):
            return (pre, -1)
    return (x,0)

class Genome(object):
    def __init__(self,fi):
        self.file = fi

    def primerMatch(self,locus,seq,ampsize):
        # get sequence with flank
        chromStart = locus.offset-ampsize[1] if locus.reverse else locus.offset+locus.length+ampsize[0]
        chromEnd   = locus.offset-ampsize[0] if locus.reverse else locus.offset+locus.length+ampsize[1]
        with pysam.FastaFile(self.file) as fasta:
            seqslice = fasta.fetch(locus.chrom,chromStart,chromEnd)
        # find sequence
        qrySeq = seq if locus.reverse else seq.translate(revcmp)[::-1]
        # create new loci
        loci = []
        for i in [ match.start() for match in re.finditer(re.escape(qrySeq), seqslice) ]:
            tm = primer3.calcTm(qrySeq)
            loci.append(Locus(locus.chrom, chromStart+i, len(qrySeq), not locus.reverse, tm))
        return loci

'''just a wrapper for pysam'''
class MultiFasta(object):
    def __init__(self,fi):
        self.file = fi
        # check sequence uniqueness
        with pysam.FastaFile(self.file) as fasta:
            self.references = fasta.references
            if len(set(fasta.references))!=len(fasta.references):
                print >> sys.stderr, self.file
                raise Exception('DuplicateSequenceNames')

    def createPrimers(self,db,bowtie='bowtie2', delete=True, tags={}, tmThreshold=50.0, endMatch=6, maxAln=20):
        # run bowtie (max 1000 alignments, allow for one gap/mismatch?)
        mapfile = self.file+'.sam'
        if not os.path.exists(mapfile):
            proc = subprocess.check_call( \
                [bowtie, '-f', '--end-to-end', '-p 2', \
                '-k '+str(maxAln), '-L 10', '-N 1', '-D 20', '-R 3', \
                '-x', db, '-U', self.file, '>', mapfile ])
        # Read fasta file (Create Primer)
        primers = {}
        with pysam.FastaFile(self.file) as fasta:
            for s in fasta.references:
                # parse target locus from fasta file
                try:
                    primername, targetposition = s.split('|')
                    reTargetposition = re.match(r'(\w+):(\d+)-(\d+):([+-])',targetposition)
                except:
                    primername = s
                    targetLocus = None
                else:
                    # create stranded targetlocus
                    reverse = True if reTargetposition.group(4)=='-' else False
                    tm = primer3.calcTm(fasta.fetch(s))  # assume targetlocus is full match
                    targetLocus = Locus(reTargetposition.group(1), int(reTargetposition.group(2)), int(reTargetposition.group(3))-int(reTargetposition.group(2)), reverse, tm)
                # create primer (with target locus)
                primertag = tags[primername] if primername in tags.keys() else None
                primers[primername] = Primer(primername,fasta.fetch(s),targetLocus,tag=primertag)

        # read SAM OUTPUT and filter alignments
        mappings = pysam.Samfile(mapfile,'r')
        alnCount = Counter()  # count alignments to kill locations of non-specific primers (count == -k)
        for aln in mappings:
            primername = aln.qname.split('|')[0]
            if aln.is_unmapped:
                continue
            else:
                alnCount[primername] += 1
            ## get reference sequence
            qry = aln.query_sequence.upper()
            ref = aln.get_reference_sequence().upper()
            refrc = ref.translate(revcmp)[::-1]
            aln_tm = primer3.calcHeterodimerTm(qry,refrc)
            # TmThreshold and mimatches in 3'end check
            if aln_tm > tmThreshold:
                if len(qry)>endMatch and len(ref)>endMatch:
                    if len([ x for x in zip(qry[-endMatch:], ref[-endMatch:]) if x[0]!=x[1] ]) == 0:
                        primers[primername].addTarget(mappings.getrname(aln.reference_id), aln.pos, aln.is_reverse, aln_tm)
        # remove primer locations for those that have hit maximum
        for k, v in primers.items():
            if len(v.loci) >= maxAln:
                v.loci = []
        # cleanup
        if delete:
            os.unlink(self.file+'.sam') # delete mapping FILE
        return primers.values()

'''Boundary exceeded exception (max list size)'''
class BoundExceedError(Exception):
    pass

'''storage location class'''
class Location(object):
    def __init__(self, vessel, well):
        # store vessel
        try:
            m = re.search(r'(\d+)',str(vessel))  # get number from vessel
            assert m
        except:
            raise
        else:
            self.vesselnumber = m.group(1)
        # store wells
        self.wells = set(well.split(','))
        try:
            assert all(map(lambda x: re.match(r'\w\d',x),list(self.wells)))
        except:
            raise Exception('InvalidWell')

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__)

    def __str__(self):
        return '-'.join([str(self.vessel()),self.well()])

    def __eq__(self,other):
        if other is None:
            return False
        return self.vesselnumber == other.vesselnumber and len(self.wells.symmetric_difference(other.wells)) == 0

    def merge(self,other):
        try:
            assert self.vesselnumber == other.vesselnumber
        except:
            raise
        else:
            self.wells.update(other.wells)

    def vessel(self):
        return int(self.vesselnumber)

    def well(self):
        return ','.join(sorted(list(self.wells)))

'''primer pair (list)'''
class PrimerPair(list):
    def __init__(self, elements, length=2, name=None, reverse=False):
        list.__init__(self, elements)
        self.length = length  # pair of primers by default
        self.reversed = reverse
        self.name = name
        self.variants = []  # list of intervals with metadata from input table
        if not name and all(self):
            commonPrefix(self[0].name, self[1].name)

    def _check_item_bound(self):
        if self.length and len(self) >= self.length:
            raise BoundExceedError()

    def _check_list_bound(self, L):
        if self.length and len(self) + len(L) > self.length:
            raise BoundExceedError()

    def reverse(self):
        super(PrimerPair, self).reverse()
        self.reversed = not self.reversed
        return self

    def append(self, x):
        self._check_item_bound()
        return super(PrimerPair, self).append(x)

    def extend(self, L):
        self._check_list_bound(L)
        return super(PrimerPair, self).extend(L)

    def insert(self, i, x):
        self._check_item_bound()
        return super(PrimerPair, self).insert(i, x)

    def __hash__(self):
        return hash(self.__str__())

    def __eq__(self,other):
        return self.name == other.name

    def __add__(self, L):
        self._check_list_bound(L)
        return super(PrimerPair, self).__add__(L)

    def __iadd__(self, L):
        self._check_list_bound(L)
        return super(PrimerPair, self).__iadd__(L)

    def __setslice__(self, *args, **kwargs):
        if len(args) > 2 and self.length:
            left, right, L = args[0], args[1], args[2]
            if right > self.length:
                if left + len(L) > self.length:
                    raise BoundExceedError()
            else:
                len_del = (right - left)
                len_add = len(L)
                if len(self) - len_del + len_add > self.length:
                    raise BoundExceedError()
        return super(PrimerPair, self).__setslice__(*args, **kwargs)

    def __lt__(self,other):
        return self.sortvalues() < other.sortvalues()

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__)

    def __str__(self):
        return '{}\t{}\t{}\t{}\t{:.1f}\t{:.1f}\t{}\t{:.1f}\t{:.1f}\t{}\t{}\t{}'.format(
            self.name,
            str(self[0].location) if self[0] and self[0].location else '',
            str(self[1].location) if self[1] and self[1].location else '',
            str(self[0].tag)+'-'+self[0].seq if self[0] else 'NO_SEQUENCE',
            self[0].tm if self[0] else 0,
            self[0].gc if self[0] else 0,
            str(self[1].tag)+'-'+self[1].seq if self[1] else 'NO_SEQUENCE',
            self[1].tm if self[1] else 0,
            self[1].gc if self[1] else 0,
            self[0].targetposition.chrom if self[0] and self[1] and self[0].targetposition else '',
            self[0].targetposition.offset+self[0].targetposition.length if self[0] and self[1] and self[0].targetposition else '',
            self[1].targetposition.offset if self[0] and self[1] and self[1].targetposition else '')

    def sequencingTarget(self):
        return (self[0].targetposition.chrom if self[0] and self[1] and self[0].targetposition else None, \
                self[0].targetposition.offset+self[0].targetposition.length if self[0] and self[1] and self[0].targetposition else None,
                self[1].targetposition.offset if self[0] and self[1] and self[1].targetposition else None,
                self.reversed)


    def targetLength(self,includePrimers=False):
        if self[0] and self[1] and self[0].targetposition and self[0].targetposition:
            if includePrimers:
                return self[1].targetposition.offset+self[1].targetposition.length-self[0].targetposition.offset
            else:
                return self[1].targetposition.offset-(self[0].targetposition.offset+self[0].targetposition.length)
        return None

    def locations(self):
        return [ self[0].location if self[0] else None, self[1].location if self[1] else None ]

    def TmRange(self):
        tm = sorted([ float(p.tm) if p.tm else 0 for p in self ])
        return '-'.join([ '{:.1f}'.format(t) for t in tm ])

    def log(self,logfile):
        timestamp = datetime.datetime.now().isoformat()
        with open(logfile, 'a') as fh:
            for p in self:
                threeprimesnps = len([ s for s in p.snp if s[1] >= 2*len(p)/3 ])
                print >> fh, '{timestamp:26} {primername:20} rank:{rank:<3d} ({fwdrev:1}) {seq:25} snps:{snps:4} misprime:{misprime:2d}'.format(\
                    timestamp=timestamp[:26], primername=p.name, rank=p.rank, fwdrev='-' if p.targetposition.reverse else '+',
                    seq=p.seq, snps=str(len(self[0].snp)-threeprimesnps)+'+'+str(threeprimesnps), misprime=len(p.loci)-1)
        return

    def pruneRanks(self):
        # set and prune ranks from name
        try:
            self[0].rank = int(self[0].name.split('_')[-2])  # if coming from primer3 (name_rank_lr)
            self[1].rank = int(self[1].name.split('_')[-2])  # if coming from primer3 (name_rank_lr)
            assert self[0].rank == self[1].rank
        except:
            raise Exception('RankError')
        else:
            ranksuffix = '_'+str(self[0].rank)
            try:
                assert self.name.endswith(ranksuffix)
            except:
                raise
            else:
                self[0].name = self.name[:-len(ranksuffix)] + self[0].name[len(self.name):]
                self[1].name = self.name[:-len(ranksuffix)] + self[1].name[len(self.name):]
        # update own name if needed
        if self.name not in self[0].name or self.name not in self[1].name:
            self.name = commonPrefix(self[0].name, self[1].name)
        return

    def sortvalues(self):
        assert len(self)==2
        return (len(self.amplicons([0,10000]))-1, self.criticalsnp(), self.mispriming(), self.snpcount(), self.designrank())

    def amplicons(self, sizeRange=[0,10000], autoreverse=True):  # counts possible amplicons to a certain size
        amplicons = []
        for m in self[0].loci:
            for n in self[1].loci:
                if m.chrom == n.chrom:
                    amplen = n.offset + n.length - m.offset
                    if (not sizeRange) or (amplen >= sizeRange[0] and amplen <= sizeRange[1]):
                        amp = (m, n, Interval(m.chrom,m.offset,n.offset + n.length,self.name))
                        amplicons.append(amp)
        if autoreverse and not amplicons:  # reverse if and find amplicons if one direction doesnt yield any
            for m in self[1].loci:
                for n in self[0].loci:
                    if m.chrom == n.chrom:
                        amplen = n.offset + n.length - m.offset
                        if (not sizeRange) or (amplen >= sizeRange[0] and amplen <= sizeRange[1]):
                            amp = (m, n, Interval(m.chrom,m.offset,n.offset + n.length,self.name))
                            amplicons.append(amp)
            if amplicons:
                self.reverse()
        return amplicons

    def snpcount(self):
        return len(self[0].snp)+len(self[1].snp)

    def mispriming(self):
        return max(max(len(self[0].loci), len(self[1].loci))-1,0)

    def criticalsnp(self):
        return len([ s for s in self[0].snp if s[1] >= 2*len(self[0])/3 ]) + \
            len([ s for s in self[1].snp if s[1]+s[2] <= len(self[1])/3 ])

    def designrank(self):
        assert self[0].rank == self[1].rank
        return int(self[0].rank)

    def check(self, limits):
        for k,v in limits.items():
            x = getattr(self,k)()
            try:
                x = int(x)
            except TypeError:
                x = len(x)
            except:
                raise
            if x > v:
                return False
        return True

    def uniqueid(self):
        return sha1(','.join([str(self[0].tag)+'-'+self[0].seq,str(self[1].tag)+'-'+self[1].seq])).hexdigest()

    '''returns primer suffixes'''
    def primerSuffixes(self):
        suffixes = []
        for p in self:
            try:
                suffix1 = '...'+p.name[len(self.name):] if p.name.startswith(self.name) else p.name
            except:
                suffixes.append('NA')
            else:
                suffixes.append(suffix1)
        return tuple(suffixes)

    '''changes name to any longer common primer name prefix'''
    def fixName(self):
        firstDifferent = min([ i for i,x in enumerate(zip(self[0].name,self[1].name)) if len(set(x))!=1 ])
        newName = self[0].name[:firstDifferent].rstrip('_-')
        if newName != self.name and len(newName) >= len(self.name):
            print >> sys.stderr, 'INFO: Renamed PrimerPair {} -> {}'.format(self.name, newName)
            self.name = newName
            return True
        return False

    def rename(self,renamer):
        oldname = self.name
        newname = renamer if type(renamer) is str else renamer(self.name)
        self.name = re.sub('^'+oldname,newname,self.name)
        for i in range(len(self)):
            self[i].name = re.sub('^'+oldname,newname,self[i].name)
        return

'''fasta/primer'''
class Primer(object):
    def __init__(self,name,seq,targetposition=None,tag=None,loci=[],location=None):
        self.rank = -1
        self.name = name
        self.seq = str(seq.upper())
        self.tag = tag
        self.tm = primer3.calcTm(self.seq)
        self.gc = (self.seq.count('G') + self.seq.count('C')) / float(len(self.seq))
        self.loci = []  # genome matches
        self.snp = []  # same order as loci attribute
        self.meta = {}  # metadata
        self.targetposition = targetposition
        self.location = location  # storage location
        if loci:
            pass

    def __hash__(self):
        return hash(self.name) ^ hash(self.seq) ^ hash(self.tag)

    def __repr__(self):
        return '<Primer ('+self.name+'):'+str(self.tag)+'-'+self.seq+' Mappings:'+str(len(self.loci))+' Target:'+str(self.targetposition)+'>'
        #return "%s(%r)" % (self.__class__, self.__dict__)

    def __str__(self):
        return '{:<20}\t{:>}-{:<}\t{:>}\t{:.2f}\t{:.1f}\t{:<}'.format(\
            self.name, str(self.tag), self.seq, str(self.location), self.tm, self.gc, str(self.targetposition))

    def __len__(self):
        return len(self.seq)

    def snpFilter(self,position):
        # return list of boolean if spliced position has Variant
        for l in self.loci:
            raise NotImplementedError

    def fasta(self,seqname=None):
        if not seqname:
            seqname = self.name
        if 'POSITION' in self.meta.keys():  # append locus
            strand = "-" if self.name.endswith('RIGHT') else '+'
            seqname += '|'+self.meta['POSITION'][0]+':'+"-".join(map(str,self.meta['POSITION'][1:]))+':'+strand
        return "\n".join([ ">"+seqname, self.seq ])

    def addTarget(self, chrom, pos, reverse, tm=None):
        self.loci.append(Locus(chrom,pos,len(self),reverse,tm))
        return

    def snpCheckPrimer(self,vcf):
        self.snp = self.targetposition.snpCheck(vcf)
        return True if self.snp else False

    def checkTarget(self):
        if self.targetposition is not None:
            for locus in self.loci:
                if locus.chrom == self.targetposition.chrom:
                    if int(locus.offset) == int(self.targetposition.offset):
                        return True
        return False

'''Locus'''
class Locus(object):
    def __init__(self,chrom,offset,length,reverse,tm):
        self.chrom = chrom
        self.offset = offset
        self.length = length
        self.reverse = reverse
        self.tm = tm

    def __str__(self):
        strand = '-' if self.reverse else '+'
        return self.chrom+":"+str(self.offset)+":"+strand

    def __lt__(self,other):
        return (self.chrom, self.offset) < (other.chrom, other.offset)

    def __eq__(self,other):
        return self.chrom == other.chrom and \
            self.offset == other.offset and \
            self.length == other.length and \
            self.reverse == other.reverse

    def __hash__(self):
        return hash((self.chrom, self.offset, self.length, self.reverse))

    def snpCheck(self,database):
        db = pysam.TabixFile(database)
        try:
            snps = db.fetch(self.chrom,self.offset,self.offset+self.length)
        except ValueError:
            snps = []
        except:
            raise
        # query database and translate to primer positions
        snp_positions = []
        for v in snps:
            f = v.split()
            snpOffset = (int(f[1])-1) - self.offset  # convert to 0-based
            snpLength = max(map(len,[ f[3] ] + f[4].split(',')))
            snp_positions.append( (f[0],snpOffset,snpLength,f[2]) )
        return snp_positions

'''primer3 wrapper class'''
class Primer3(object):
    def __init__(self,genome,target,flank=200):
        self.genome = genome
        self.target = target
        self.flank = flank
        fasta = pysam.FastaFile(self.genome)
        self.designregion = ( str(self.target[0]), self.target[1]-self.flank, self.target[2]+self.flank )
        self.sequence = fasta.fetch(*self.designregion)
        self.pairs = []
        self.explain = []

    def __len__(self):
        return len(self.pairs)

    def design(self,name,pars):
        # Sequence args
        seq = {
            'SEQUENCE_ID': str(name),
            'SEQUENCE_TEMPLATE': str(self.sequence),
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [0, self.flank, len(self.sequence)-self.flank, self.flank]
        }
        # design primers
        primers = primer3.bindings.designPrimers(seq,pars)
        # parse primer
        primerdata, explain = defaultdict(dict), []
        for k,v in primers.items():
            m = re.match(r'PRIMER_(RIGHT|LEFT)_(\d+)(.*)',k)
            if m:
                primername = name+"_"+str(m.group(2))+'_'+m.group(1)
                if m.group(3):
                    primerdata[primername][m.group(3)[1:]] = v
                else:
                    absoluteStart = self.designregion[1]+v[0]-(v[1]-1) if m.group(1)=="RIGHT" else self.designregion[1]+v[0]
                    absoluteEnd = self.designregion[1]+v[0] if  m.group(1)=="RIGHT" else self.designregion[1]+v[0]+v[1]
                    primerdata[primername]['POSITION'] = (self.designregion[0], absoluteStart, absoluteEnd)
            elif k.endswith('EXPLAIN'):
                self.explain.append(v)
        designedPrimers, designedPairs = {}, {}
        for k,v in sorted(primerdata.items()):
            # k primername # v dict of metadata
            if 'SEQUENCE' not in designedPrimers.keys():
                designedPrimers[v['SEQUENCE']] = Primer(k,v['SEQUENCE'])
                m = re.search(r'(\d+)_(LEFT|RIGHT)',k)
                # store pairs (reference primers)
                if int(m.group(1)) not in designedPairs.keys():
                    designedPairs[int(m.group(1))] = PrimerPair([None, None])
                designedPairs[int(m.group(1))][0 if m.group(2).startswith('LEFT') else 1] = designedPrimers[v['SEQUENCE']]
                # all other datafields
                designedPrimers[v['SEQUENCE']].meta = v
        # store
        self.pairs = OrderedDict(sorted(designedPairs.items())).values()
        return len(self.pairs)


if __name__=="__main__":
    mf = MultiFasta(sys.argv[1])
    primers = mf.createPrimers('/Users/dbrawand/dev/snappy/WORK/genome/human_g1k_v37.bowtie')
    for primer in primers:
        print primer
