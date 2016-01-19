#!/usr/bin/env python
import sys, os, re, datetime
from hashlib import md5, sha1
import primer3
import pysam
import subprocess
from collections import defaultdict, OrderedDict
from .interval import Interval

'''returns common prefix (substring)'''
def commonPrefix(left,right,stripchars='-_ ',commonlength=3):
    matchingPositions = [ i+1 for i,j in enumerate([ i for i, x in enumerate(zip(left,right)) if len(set(x)) == 1]) if i==j]
    if matchingPositions and max(matchingPositions) >= commonlength:
        return left[:max(matchingPositions)].rstrip(stripchars)
    else:
        return None

'''just a wrapper for pysam'''
class MultiFasta(object):
    def __init__(self,file):
        self.file = file
        # check sequence uniqueness
        with pysam.FastaFile(self.file) as fasta:
            if len(set(fasta.references))!=len(fasta.references):
                print >> sys.stderr, self.file
                raise Exception('DuplicateSequenceNames')

    def createPrimers(self,db,bowtie='bowtie2'):
        # run bowtie (max 1000 alignments, allow for one gap/mismatch?)
        mapfile = self.file+'.sam'
        if not os.path.exists(mapfile):
            proc = subprocess.check_call( \
                [bowtie, '-f', '--end-to-end', \
                '-k 50', '-L 10', '-N 1', '-D 20', '-R 3', \
                '-x', db, '-U', self.file, '>', mapfile ])
        # Read fasta file (Create Primer)
        primers = {}
        fasta = pysam.FastaFile(self.file)
        for s in fasta.references:
            # parse target locus from fasta file
            try:
                primername, targetposition = s.split('|')
                reTargetposition = re.match(r'(\w+):(\d+)-(\d+)',targetposition)
            except:
                #raise Exception('PrimerNameError')
                primername = s
                targetLocus = None
            else:
                reverse = True if primername.split('_')[-1].startswith("r") else False
                targetLocus = Locus(reTargetposition.group(1), int(reTargetposition.group(2)), int(reTargetposition.group(3))-int(reTargetposition.group(2)), reverse)
            # create primer (with target locus)
            primers[primername] = Primer(primername,fasta.fetch(s),targetLocus)

        # read SAM OUTPUT
        mappings = pysam.Samfile(mapfile,'r')
        for aln in mappings:
            primername = aln.qname.split('|')[0]
            # add full matching loci
            if not any(zip(*aln.cigar)[0]): # all matches (full length)
                primers[primername].addTarget(mappings.getrname(aln.reference_id), aln.pos, aln.is_reverse)
            # add other significant matches (1 mismatch/gap)
            elif zip(*aln.cigar)[0].count(0) >= len(aln.seq)-1:
                primers[primername].sigmatch += 1

        ##os.unlink(self.file+'.sam') # delete mapping FILE
        return primers.values()


'''Boundary exceeded exception (max list size)'''
class BoundExceedError(Exception):
    pass


'''primer pair (list)'''
class PrimerPair(list):
    def __init__(self,elements,length=2, location=[]):
        list.__init__(self, elements)
        self.length = length  # pair of primers by default
        self.reversed = False

    def _check_item_bound(self):
        if self.length and len(self) >= self.length:
            raise BoundExceedError()

    def _check_list_bound(self, L):
        if self.length and len(self) + len(L) > self.length:
            raise BoundExceedError()

    def reverse(self):
        super(PrimerPair, self).reverse()
        self.reversed =  not self.reversed
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
        if self[0].targetposition and self[1].targetposition:
            return '{}\t{}\t{:.1f}\t{:.1f}\t{}\t{:.1f}\t{:.1f}\t{}\t{}\t{}'.format(self.name(), \
                self[0].seq, self[0].tm, self[0].gc, \
                self[1].seq, self[1].tm, self[1].gc, \
                self[0].targetposition.chrom, self[0].targetposition.offset+self[0].targetposition.length, self[1].targetposition.offset)
        else:
            return '{}\t{}\t{:.1f}\t{:.1f}\t{}\t{:.1f}\t{:.1f}\tNA\tNA\tNA'.format(self.name(), \
                self[0].seq, self[0].tm, self[0].gc, \
                self[1].seq, self[1].tm, self[1].gc)

    def log(self,logfile):
        with open(logfile, 'a') as fh:
            print >> fh, "PRIMER PAIR:", self.name()
            # date, primer name, FWD/REV,
            # sequence, 
            # snpcount, mispriming, criticalsnp, designrank
            if self[0].targetposition.reverse:
                leftReverse = 'REV'
            else:
                leftReverse = 'FWD'

            if self[1].targetposition.reverse:
                rightReverse = 'REV'
            else:
                rightReverse = 'FWD'

            print >> fh, '\n\t', "Design time:", datetime.datetime.now(), " |  Primer name:", self[0].name, \
            " |  Fwd/Rev:", leftReverse
            print >> fh, '\n\t', "Sequence:", self[0].seq


            print >> fh, '\n\n\t', "Design time:", datetime.datetime.now(), " |  Primer name:", self[1].name, \
            " |  Fwd/Rev:", rightReverse
            print >> fh, '\n\t', "Sequence:", self[1].seq

            print >> fh, '\n\t', "SNPcount:", self.snpcount(), " |  Misprimes:", self.mispriming(), \
            " |  Critical SNPs:", self.criticalsnp(), " |  Primer3 design rank:", self.designrank(), '\n'
           
            print >> fh, '-'*41


        return


    def name(self):
        return commonPrefix(self[0].name, self[1].name)

    def pruneRanks(self):
        # set and prune ranks from name
        try:
            self[0].rank = int(self[0].name.split('_')[-2])  # if coming from primer3 (name_rank_lr)
            self[1].rank = int(self[1].name.split('_')[-2])  # if coming from primer3 (name_rank_lr)
            assert self[0].rank == self[1].rank
        except:
            raise Exception('RankError')
        else:
            pairname = self.name()
            ranksuffix = '_'+str(self[0].rank)
            try:
                assert pairname.endswith(ranksuffix)
            except:
                raise
            else:
                self[0].name = pairname[:-len(ranksuffix)] + self[0].name[len(pairname):]
                self[1].name = pairname[:-len(ranksuffix)] + self[1].name[len(pairname):]
        return

    def sortvalues(self):
        assert len(self)==2
        return (len(self.amplicons([0,10000]))-1, self.criticalsnp(), self.mispriming(), self.snpcount(), self.designrank())

    def amplicons(self,sizeRange=[0,1000]):  # counts possible amplicons to a certain size
        amplicons = []
        for m in self[0].loci:
            for n in self[1].loci:
                if m.chrom == n.chrom:
                    amplen = n.offset + n.length - m.offset
                    if (not sizeRange) or (amplen >= sizeRange[0] and amplen <= sizeRange[1]):
                        amp = (m, n, Interval(m.chrom,m.offset,n.offset + n.length,self.name()))
                        amplicons.append(amp)
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
        return sha1(','.join([self[0].seq,self[1].seq])).hexdigest()

    def __hash__(self):
        return hash(self[0]) ^ hash(self[1])


'''fasta/primer'''
class Primer(object):
    def __init__(self,name,seq,targetposition=None,tm=None,gc=None,loci=[]):
        self.rank = -1
        self.name = name
        self.seq = str(seq.upper())
        self.tm = tm
        self.gc = gc
        self.loci = []  # genome matches
        self.snp = []  # same order as loci attribute
        self.meta = {}  # metadata
        self.sigmatch = 0  # significant other matches (just counted)
        self.targetposition = targetposition
        if loci:
            pass

    def __str__(self):
        return '<Primer ('+self.name+'): '+self.seq+', Targets: '+str(len(self.loci))+' (other significant: '+str(self.sigmatch)+'); Target position: '+str(self.targetposition)+'>'

    def __repr__(self):
        return '{:<20} {:<24} {:<}'.format(self.name,self.seq,str(self.targetposition))

    def __len__(self):
        return len(self.seq)

    def snpFilter(self,position):
        # return list of boolean if spliced position has Variant
        for l in self.loci:
            raise NotImplementedError

    def fasta(self,seqname=None):
        if not seqname:
            seqname = self.name
        if 'POSITION' in self.meta.keys():
            seqname += '|'+self.meta['POSITION'][0]+':'+"-".join(map(str,self.meta['POSITION'][1:]))
        return "\n".join([ ">"+seqname, self.seq ])

    def addTarget(self, chrom, pos, reverse):
        self.loci.append(Locus(chrom,pos,len(self),reverse))
        return

    def calcProperties(self):
        # get Tm via primer3
        self.tm = primer3.calcTm(self.seq)
        # calc GC
        self.gc = (self.seq.count('G') + self.seq.count('C')) / float(len(self.seq))
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
    def __init__(self,chrom,offset,length,reverse):
        self.chrom = chrom
        self.offset = offset
        self.length = length
        self.reverse = reverse

    def __str__(self):
        strand = '-' if self.reverse else '+'
        return self.chrom+":"+str(self.offset)+":"+strand

    def __lt__(self,other):
        return (self.chrom, self.offset) < (other.chrom, other.offset)

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
        self.designregion = ( self.target[0], self.target[1]-self.flank, self.target[2]+self.flank )
        self.sequence = fasta.fetch(*self.designregion)
        self.pairs = []
        self.explain = []

    def __len__(self):
        return len(self.pairs)

    def design(self,name,pars):
        # extract sequence with flanks
        # Sequence args
        seq = {
            'SEQUENCE_ID': name,
            'SEQUENCE_TEMPLATE': self.sequence,
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [0, self.flank, len(self.sequence)-self.flank, self.flank]
        }
        # design primers
        primers = primer3.bindings.designPrimers(seq,pars)
        # parse primer
        primerdata, explain = defaultdict(dict), []
        for k,v in primers.items():
            # print k, v
            # print self. designregion
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
                designedPrimers[v['SEQUENCE']].calcProperties()

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
