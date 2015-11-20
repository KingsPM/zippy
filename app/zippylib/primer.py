#!/usr/bin/env python
import sys, os, re
from hashlib import md5
import primer3
import pysam
import subprocess
from collections import defaultdict, OrderedDict


'''just a wrapper for pysam'''
class MultiFasta(object):
    def __init__(self,file):
        self.file = file

    def createPrimers(self,db,bowtie='bowtie2'):
        # run bowtie (max 1000 alignments, allow for one gap/mismatch?)
        mapfile = self.file+'.sam'
        if not os.path.exists(mapfile):
            proc = subprocess.check_call( \
                [bowtie, '-f', '--end-to-end', \
                '-k 10', '-L 10', '-N 1', '-D 20', '-R 3', \
                '-x', db, '-U', self.file, '>', mapfile ])
        # read SAM OUTPUT
        primers = {}
        mappings = pysam.Samfile(mapfile,'r')
        print self.file
        for aln in mappings:
            #print aln.rname, aln.qname, aln.pos, aln.seq
            if aln.qname not in primers.keys():
                # create primer
                primers[aln.qname] = Primer(aln.qname,aln.seq)
            # add full matching loci
            if not any(zip(*aln.cigar)[0]): # all matches (full length)
                primers[aln.qname].addTarget(mappings.getrname(aln.reference_id), aln.pos, aln.is_reverse)
            # add other significant matches (1 mismatch/gap)
            elif zip(*aln.cigar)[0].count(0) >= len(aln.seq)-1:
                primers[aln.qname].sigmatch += 1

        ## delete mapping FILE
        ####os.unlink(self.file+'.sam')
        return primers.values()

'''fasta/primer'''
class Primer(object):
    def __init__(self,name,seq,tm=None,gc=None,loci=[]):
        self.name = name if name else 'primer_'+md5(seq).hexdigest()[:8]
        self.seq = seq.upper()
        self.tm = tm
        self.gc = gc
        self.loci = []  # genome matches
        self.snp = []  # same order as loci attribute
        self.meta = {}  # metadata
        self.sigmatch = 0  # significant other matches (just counted)
        if loci:
            pass

    def __str__(self):
        return '<Primer ('+self.name+'): '+self.seq+', Targets: '+str(len(self.loci))+' (other significant: '+str(self.sigmatch)+')>'

    def __repr__(self):
        return '<Primer ('+self.name+'): '+self.seq+', Targets: '+str(len(self.loci))+' (other significant: '+str(self.sigmatch)+')>'

    def __len__(self):
        return len(self.seq)

    # def __repr__(self):
    #     '''primer sequence with locations and annotations'''
    #     raise NotImplementedError

    def snpFilter(self,position):
        # return list of boolean if spliced position has Variant
        for l in self.loci:
            raise NotImplementedError


    def fasta(self,seqname=None):
        if not seqname:
            seqname = self.name
        if 'POSITION' in self.meta.keys():
            seqname += ' '+self.meta['POSITION'][0]+':'+"-".join(map(str,self.meta['POSITION'][1:]))
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


'''Locus'''
class Locus(object):
    def __init__(self,chrom,offset,length,reverse):
        self.chrom = chrom
        self.offset = offset
        self.length = length
        self.reverse = reverse
        self.snps = []

    def __str__(self):
        strand = '-' if self.reverse else '+'
        return self.chrom+":"+str(self.offset)+":"+strand

    def __lt__(self,other):
        return (self.chrom, self.offset) < (other.chrom, other.offset)

    def snpCheck(self,database):
        db = pysam.TabixFile(database)
        print '\tLOCUS', str(self), db
        try:
            snps = db.fetch(self.chrom,self.offset,self.offset+self.length)
        except ValueError:
            snps = []
        except:
            raise
        # query database and translate to primer positions
        for v in snps:
            print "\t\tSNP", v 
            f = v.split()
            snpOffset = int(f[1])-self.offset
            snpLength = max(len(f[3]),len(f[4]))
            self.snps.append( (f[0],snpOffset,snpLength,f[2]) )
        return


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
            m = re.match(r'PRIMER_(RIGHT|LEFT)_(\d+)(.*)',k)
            if m:
                primername = name+"_"+str(m.group(2))+'_'+m.group(1)
                if m.group(3):
                    primerdata[primername][m.group(3)[1:]] = v
                else:
                    primerdata[primername]['POSITION'] = (self.designregion[0], self.designregion[1]+v[0], self.designregion[1]+v[0]+v[1])
            elif k.endswith('EXPLAIN'):
                self.explain.append(v)

        designedPrimers, designedPairs = {}, {}
        for k,v in sorted(primerdata.items()):
            # k primername
            # v dict of metadata 
            if 'SEQUENCE' not in designedPrimers.keys():
                designedPrimers[v['SEQUENCE']] = Primer(k,v['SEQUENCE'])  # no name autosets
                designedPrimers[v['SEQUENCE']].calcProperties()

                m = re.search(r'(\d+)_(LEFT|RIGHT)',k)
                # store pairs (reference primers)
                if int(m.group(1)) not in designedPairs.keys():
                    designedPairs[int(m.group(1))] = [None, None]
                designedPairs[int(m.group(1))][0 if m.group(2).startswith('LEFT') else 1] = designedPrimers[v['SEQUENCE']]
                # all other datafields
                designedPrimers[v['SEQUENCE']].meta = v
        # store
        self.pairs += OrderedDict(sorted(designedPairs.items())).values()
        # if design fails there will be 0 pairs, simple!
        # if not self.pairs:
        #     print self.explain
        #     print self.sequence
        #     #raise Exception('DesignFail')
        return len(self.pairs)

    def show(self,width=100):
        # get features
        for p in self.pairs:
            f = { self.flank:'[', len(self.sequence)-self.flank:']' }
            f[int(p[0].meta[''][0] + p[0].meta[''][1])]= '>'
            f[int(p[1].meta[''][0])] = '<'
            # get fstring positions
            fpos = defaultdict(list)
            for k in sorted(f.keys()):
                p = int(k*width/float(2*self.flank+(self.target[2]-self.target[1])))
                fpos[p].append(f[k])
            # create featurestring
            fstring = ''
            for p in sorted(fpos.keys()):
                # spacer
                if len(fstring) < p:
                    fstring += ' ' * (p-len(fstring))
                # char
                if len(fpos[p]) > 1:
                    fstring += "*"
                else:
                    fstring += fpos[p][0]
            # fill end
            if len(fstring) < width:
                fstring += ' ' * (width-len(fstring))
            # print
            print self.target[0], self.target[1],'|', fstring,'|', self.target[2]
        return

if __name__=="__main__":
    mf = MultiFasta(sys.argv[1])
    primers = mf.createPrimers('/Users/dbrawand/dev/snappy/WORK/genome/human_g1k_v37.bowtie')
    for primer in primers:
        print primer
