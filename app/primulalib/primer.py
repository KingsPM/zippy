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
        ## run bowtie (max 100 alignments, allow for one gap/mismatch?)
        ####proc = subprocess.check_call( \
        ####    [bowtie, '-f', '--end-to-end', \
        ####    '-k 100', '-L 10', '-N 1', '-D 20', '-R 3', \
        ####    '-x', db, '-U', self.file, '>', self.file+'.sam' ])
        # read SAM OUTPUT
        primers = {}
        mappings = pysam.Samfile(self.file+'.sam','r')
        for aln in mappings:
            #print aln.rname, aln.qname, aln.pos, aln.seq
            if aln.qname not in primers.keys():
                # create primer
                primers[aln.qname] = Primer(aln.qname,aln.seq)
            # add full matching loci
            if not any(zip(*aln.cigar)[0]): # all matches (full length)
                primers[aln.qname].addTarget(aln.rname, aln.pos, aln.is_reverse)
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
        self.loci = []
        self.meta = {}  # metadata
        self.sigmatch = 0  # significant other matches (just counted)
        if loci:
            pass

    def __str__(self):
        return '<Primer ('+self.name+'): '+self.seq+', Targets: '+str(len(self.loci))+' (other significant: '+str(self.sigmatch)+')>'

    def addTarget(self, chrom, pos, reverse):
        self.loci.append([chrom,pos,reverse])
        return

    def calcProperties(self):
        # get Tm via primer3
        self.tm = primer3.calcTm(self.seq)
        # calc GC
        self.gc = (self.seq.count('G') + self.seq.count('C')) / float(len(self.seq))
        return


class Primer3(object):
    def __init__(self,genome,target,flank=200):
        self.genome = genome
        self.target = target
        self.flank = flank
        fasta = pysam.FastaFile(self.genome)
        flanked = ( self.target[0], self.target[1]-self.flank, self.target[2]+self.flank )
        self.sequence = fasta.fetch(*flanked)
        self.pairs = []

    def design(self,name,pars):
        # extract sequence with flanks
        # Sequence args
        seq = {
            'SEQUENCE_ID': name,
            'SEQUENCE_TEMPLATE': self.sequence,
            'SEQUENCE_INCLUDED_REGION': [0,len(self.sequence)],
            'SEQUENCE_EXCLUDED_REGION': [self.flank, len(self.sequence)-self.flank]
        }
        # design primers
        primers = primer3.bindings.designPrimers(seq,pars)
        # create primers
        designedPrimers, designedPairs = {}, {}
        for k,v in primers.items():
            m = re.match(r'PRIMER_(RIGHT|LEFT)_(\d+)_SEQUENCE',k)
            if m:
                # create primer
                if v not in designedPrimers.keys():
                    designedPrimers[v] = Primer(name+"_"+str(m.group(2))+m.group(1),v)  # no name autosets
                    designedPrimers[v].calcProperties()
                # store pairs (reference primers)
                if int(m.group(2)) not in designedPairs.keys():
                    designedPairs[int(m.group(2))] = [None, None]
                designedPairs[int(m.group(2))][0 if m.group(1).startswith('LEFT') else 1] = designedPrimers[v]
        # add metadata
        for k,v in primers.items():
            m = re.match(r'PRIMER_(RIGHT|LEFT)_(\d+)',k)
            if m:
                metavar = '_'.join(k.split('_')[3:])
                designedPairs[int(m.group(2))][0 if m.group(1).startswith('LEFT') else 1].meta[metavar] = v
        # store
        self.pairs += OrderedDict(sorted(designedPairs.items())).values()
        return len(self.pairs)

    def show(self,fh=sys.stderr,width=50):
        # get features
        for p in self.pairs:
            f = { self.flank:'[', len(self.sequence)-self.flank:']' }
            f[int(p[0].meta[''][0] + p[0].meta[''][1])]= '>'
            f[int(p[1].meta[''][0])] = '<'
            # create featurestring
            fstring = ''
            for k in sorted(f.keys()):
                p = int(k*width/float(self.target[2]-self.target[1]))
                assert p >= len(fstring)
                # spacer
                if p == len(fstring):
                    fstring += "*"
                else:
                    fstring += ' ' * (p-len(fstring))
                # char
                fstring += f[k]
            if len(fstring) < width:
                fstring += ' ' * len(fstring)-width
            # print
            print >> fh, ('{:<5} {:>10d} {:'+str(width)+'} {:<100d}').format(
                self.target[0],self.target[1],fstring,self.target[2])
        return

if __name__=="__main__":
    mf = MultiFasta(sys.argv[1])
    primers = mf.createPrimers('/Users/dbrawand/dev/snappy/WORK/genome/human_g1k_v37.bowtie')
    for primer in primers:
        print primer
