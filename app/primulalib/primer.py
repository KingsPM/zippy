#!/usr/bin/env python
import sys, os
import primer3
import pysam
import subprocess

'''just a wrapper for pysam'''
class MultiFasta(object):
    def __init__(self,file):
        self.file = file

    def createPrimers(self,db,bowtie='bowtie2'):
        ## run bowtie
        ####proc = subprocess.check_call( \
        ####    [bowtie, '-f', '--end-to-end', '-x', db, '-U', self.file, '>', self.file+'.sam' ])
        ## read SAM OUTPUT
        primers = {}
        mappings = pysam.Samfile(self.file+'.sam','r')
        for aln in mappings:
            print aln.rname, aln.qname, aln.pos, aln.seq

            if aln.qname not in primers.keys():
                # create primer
                primers[aln.qname] = Primer(aln.qname,aln.seq)
            # add full matching loci
            if not any(zip(*aln.cigar)[0]): # all matches (full length)
                primers[aln.qname].addTarget(aln.rname, aln.pos, aln.is_reverse)

        ## delete mapping FILE
        ######os.unlink(self.file+'.sam')
        return primers.values()

'''fasta/primer'''
class Primer(object):
    def __init__(self,name,seq,tm=None,gc=None,loci=[]):
        self.name = name
        self.seq = seq.upper()
        self.tm = tm
        self.gc = gc
        self.loci = []
        if loci:
            pass

    def __str__(self):
        return '<Primer: '+self.seq+', Targets: '+str(len(self.loci))+'>'

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
    def __init__(self,genome,seqflank=1000,varflank=10):
        self.genome = genome
        self.seqflank = seqflank
        self.varflank = varflank

    def design(self,name,locus,pars):
        # extract sequence with flanks
        fasta = pysam.FastaFile(self.genome)
        flanked = ( locus[0], locus[1]-self.seqflank, locus[2]+self.seqflank )
        sequence = fasta.fetch(*flanked)
        print sequence
        # Sequence args
        seq = {
            'SEQUENCE_ID': name,
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_INCLUDED_REGION': [self.seqflank-self.varflank,len(sequence)-self.seqflank+self.varflank]
        }
        # design primers
        primers = primer3.bindings.designPrimers(seq,pars)
        # parse primers

        ....

        return primers

if __name__=="__main__":
    mf = MultiFasta(sys.argv[1])
    primers = mf.createPrimers('/Users/dbrawand/dev/snappy/WORK/genome/human_g1k_v37.bowtie')
    for primer in primers:
        print primer
