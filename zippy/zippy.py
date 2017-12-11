#!/usr/bin/env python

__doc__=="""
################################################################
# Zippy - Primer database and automated design                 #
# -- Organisation: Viapath Analytics / King's College Hospital #
# -- From: 26/08/2015                                          #
################################################################
"""
__author__ = "David Brawand"
__credits__ = ['David Brawand','Christopher Wall']
__license__ = "MIT"
__version__ = "2.3.4"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import os
import re
import sys
import json
import tempfile
import hashlib
import csv
from zippylib.files import VCF, BED, GenePred, Interval, Data, readTargets, readBatch
from zippylib.primer import Genome, MultiFasta, Primer3, Primer, PrimerPair, Location, parsePrimerName
from zippylib.reports import Test
from zippylib.database import PrimerDB
from zippylib.interval import IntervalList
from zippylib import ConfigError, Progressbar, banner
from zippylib.reports import Worksheet
from argparse import ArgumentParser
from copy import deepcopy
from collections import defaultdict, Counter
from urllib import unquote
import cPickle as pickle

'''file MD5'''
def fileMD5(fi, block_size=2**20):
    md5 = hashlib.md5()
    with open(fi,'rb') as fh:
        while True:
            data = fh.read(block_size)
            if not data:
                break
            md5.update(data)
    return md5.hexdigest()

'''import primer locations from table'''
def importPrimerLocations(inputfile):
    primerlocs = {}
    with open(inputfile) as infh:
        for i,line in enumerate(infh):
            if i == 0:
                header = map(lambda x : x.lower(), line.rstrip().split('\t'))
            else:
                f = map(lambda x: x.strip('"'), line.rstrip().split('\t'))
                l = dict(zip(header,f))
                # store metadata and write fasta
                if 'vessel' in l.keys() and 'well' in l.keys() and \
                    l['vessel'] and l['well'] and l['vessel'] != 'None' and l['well'] != 'None':
                    # store location
                    loc = Location(l['vessel'],l['well'])
                    if l['primername'] in primerlocs.keys():
                        primerlocs[l['primername']].merge(loc)
                    else:
                        primerlocs[l['primername']] = loc
    return primerlocs

'''converts comma seperated variant descriptor to pseudo human readable format'''
def shortHumanReadable(x):
    fields = unquote(x).split(',')
    return '_'.join([ fields[0], hashlib.sha1(','.join(fields[:-1])).hexdigest()[:6].upper() ])

'''reads fasta/tab inputfile, searches genome for targets, decides if valid pairs'''
def importPrimerPairs(inputfile,config,primer3=True):
    # read table/fasta
    primersets = defaultdict(list)  # pair primersets
    primertags = {}  # primer tags from table
    if not inputfile.split('.')[-1].startswith('fa'):  # ignores duplicate sequence
        primerseqs = {}
        fastafile = 'import_' + fileMD5(inputfile)[:8] + '.fasta'
        with open(fastafile,'w') as outfh:
            with open(inputfile) as infh:
                for i,line in enumerate(infh):
                    if i == 0:
                        minimalHeader = set(['primername','primerset','tag','sequence','vessel','well'])
                        header = map(lambda x : x.lower(), line.rstrip().split('\t'))
                        try:
                            assert not minimalHeader.difference(set(header))
                        except:
                            print >> sys.stderr, 'ERROR: Missing columns (%s)' % ','.join(list(minimalHeader.difference(set(header))))
                            raise Exception('FileHeaderError')
                    else:
                        f = map(lambda x: x.strip('"'), line.rstrip().split('\t'))
                        l = dict(zip(header,f))
                        # remove tag from sequence
                        if l['tag']:
                            try:
                                tagseqs = config['ordersheet']['sequencetags'][l['tag']]['tags']
                            except:
                                pass
                            else:
                                for t in tagseqs:
                                    if l['sequence'].startswith(t):
                                        l['sequence'] = l['sequence'][len(t):]
                                        break
                        # store metadata and write fasta
                        if l['primername'] in primerseqs.keys():
                            try:
                                assert l['sequence'] == primerseqs[l['primername']]
                                assert l['tag'] == primertags[l['primername']]
                            except:
                                print >> sys.stderr, l['primername']
                                print >> sys.stderr, primerseqs[l['primername']]
                                print >> sys.stderr, primertags[l['primername']]
                                raise Exception('ImportFormattingError')
                        else:
                            print >> outfh, '>'+l['primername']
                            print >> outfh, l['sequence']
                        if l['primerset']:
                            primersets[l['primername']].append(l['primerset'])
                        primertags[l['primername']] = l['tag']
                        primerseqs[l['primername']] = l['sequence']
        primerfile = MultiFasta(fastafile)
    else:
        primerfile = MultiFasta(inputfile)
        # set default tags for import
        for r in primerfile.references:
            primertags[r] = config['import']['tag']
    print >> sys.stderr, "Placing primers on genome..."
    # Align primers to genome
    primers = primerfile.createPrimers(config['design']['bowtieindex'], \
        delete=False,tags=primertags, \
        tmThreshold=config['design']['mispriming']['minimaltm'], \
        endMatch=config['design']['mispriming']['identity3prime'])  # places in genome
    # pair primers (by name or by primerset) MAKE COPIES!!!!
    pairs = {}
    for p in primers:
        setnames = primersets[p.name] \
            if p.name in primersets.keys() and len(primersets[p.name])>0 \
            else [ parsePrimerName(p.name)[0] ]
        for setname in setnames:
            try:
                pairs[setname]
            except KeyError:
                try:
                    pairs[setname] = PrimerPair([None,None],name=setname)
                except:
                    print >> sys.stderr, '>>',primersets[p.name], '|', p.name, '|', setnames, '<'
                    raise
            except:
                raise
            # get primer orientation (might be wrong if guesses from name, will correct after)
            ## this basically just makes sure primers get paired (one fwd, one reverse)
            reverse = p.targetposition.reverse if p.targetposition else parsePrimerName(p.name)[1] < 0
            try:
                if reverse and pairs[setname][1] is None:
                    pairs[setname][1] = deepcopy(p)
                else:
                    if pairs[setname][0] is None:
                        pairs[setname][0] = deepcopy(p)
                    else:
                        assert pairs[setname][1] is None
                        pairs[setname][1] = deepcopy(p)
            except:
                print >> sys.stderr, "ERROR: Primer pair strand conflict?"
                print >> sys.stderr, "PRIMER0", pairs[setname][0]
                print >> sys.stderr, "PRIMER1", pairs[setname][1]
                print >> sys.stderr, "REVERSE", reverse
                print >> sys.stderr, "SETNAME", setname
                print >> sys.stderr, "PRIMER", p.name, parsePrimerName(p.name)
                print >> sys.stderr, "PAIRS", pairs[setname]
                raise
    # check if any unpaired primers
    for k,v in pairs.items():
        if not all(v):
            print >> sys.stderr, "WARNING: primer set %s is incomplete and skipped" % k
            del pairs[k]
    # prune ranks in primer3 mode (renames pair)
    if primer3:
        for p in pairs.values():
            assert p[0].targetposition and p[1].targetposition  # make sure target postiions are set
            p.pruneRanks()
        validPairs = pairs.values()
    else:  # guess target if not set
        validPairs = []
        print >> sys.stderr, 'Identifying correct amplicons for unplaced primer pairs...'
        for p in pairs.values():
            if not p[0].targetposition or not p[1].targetposition:
                amplicons = p.amplicons(config['import']['ampliconsize'],autoreverse=True)
                if amplicons:
                    shortest = sorted(amplicons,key=lambda x: len(x[2]))[0]  # sort amplicons by size
                    if len(amplicons)>1:
                        print >> sys.stderr, 'WARNING: multiple amplicons for {}. Assuming shortest ({}bp) is correct.'.format(p.name,str(len(shortest[2])))
                    p[0].targetposition = shortest[0]  # m
                    p[1].targetposition = shortest[1]  # n
                    validPairs.append(p)
                elif not primer3:
                    # try to find amplicon by sequence matching if no amplicons from genome mapping with sufficient Tm
                    refGenome = Genome(config['design']['genome'])
                    # get new loci (one round)
                    newLoci = [[],[]]
                    for mapped,query in [[0,1],[1,0]]:
                        for l in p[mapped].loci:
                            newLoci[query] += refGenome.primerMatch(l,p[query].seq,config['import']['ampliconsize'])
                    # add new loci
                    if not newLoci[0] and not newLoci[1]:
                        print >> sys.stderr, 'WARNING: {} is not specific and not imported ({},{})'.format(p.name,len(p[0].loci),len(p[1].loci))
                        continue
                    else:  # add new loci
                        for i,loc in enumerate(newLoci):
                            p[i].loci += loc
                            p[i].loci = list(set(p[i].loci))  # remove redundancy
                    # store new amplicon
                    amplicons = p.amplicons(config['import']['ampliconsize'],autoreverse=True)
                    if amplicons:
                        p[0].targetposition = amplicons[0][0]  # m
                        p[1].targetposition = amplicons[0][1]  # n
                        validPairs.append(p)
                    else:
                        print >> sys.stderr, '\n'.join([ "-> "+str(l) for l in p[0].loci])
                        print >> sys.stderr, '\n'.join([ "<- "+str(l) for l in p[1].loci])
                        print >> sys.stderr, 'WARNING: Primer set {} has no valid amplicons ({},{})'.format(p.name,len(p[0].loci),len(p[1].loci))
                else:
                    print >> sys.stderr, 'WARNING: Primer set {} does not produce a well-sized, unique amplicon ({},{})'.format(p.name,len(p[0].loci),len(p[1].loci))
            else:
                validPairs.append(p)
    return validPairs

'''get primers from intervals'''
def getPrimers(intervals, db, design, config, tiers=[0], rename=None, compatible=False):
    ivpairs = defaultdict(list)  # found/designed primer pairs (from database or design)
    blacklist = db.blacklist() if db else []
    try:
        blacklist += pickle.load(open(config['blacklistcache'],'rb'))
    except:
        print >> sys.stderr, 'Could not read blacklist cache, check permissions'
        print >> sys.stderr, os.getcwd(), config['blacklistcache']
    seqhash = lambda x,y: hashlib.sha1(','.join([x,y])).hexdigest()  # sequence pair hashing function
    # build gap primers and hash valid pairs
    if compatible:
        assert len(intervals)==2
        # get compatible pairs (no SNPcheck, just hash seqs)
        compatible = set()
        # build gap sequence
        gap = Interval(intervals[0].chrom, intervals[0].chromEnd, intervals[1].chromStart)
        # design primers
        progress = Progressbar(len(tiers),'Building compatibility list')
        for tier in tiers:
            sys.stderr.write('\r'+progress.show(tier))
            # get sequence
            maxFlank = max([ max(x) for x in config['design']['primer3'][tier]['PRIMER_PRODUCT_SIZE_RANGE'] ])
            p3 = Primer3(config['design']['genome'], gap.locus(), maxFlank)
            # alter configuration for gap pcr
            cfg = deepcopy(config['design']['primer3'][tier])
            # redefine PRIMER_PRODUCT_SIZE_RANGE
            cfg['PRIMER_PRODUCT_SIZE_RANGE'] = [ [x[0]+len(gap), x[1]+len(gap) ] for x in cfg['PRIMER_PRODUCT_SIZE_RANGE']]
            # redefine primer count (Get 2x more for compatibiliy pairs as they are cheap to design)
            cfg['PRIMER_NUM_RETURN'] *= 2
            p3.design('gap-PCR', cfg)
            if p3.pairs:
                for p in p3.pairs:
                    compatible.add(seqhash(p[0].seq, p[1].seq))
        sys.stderr.write('\r'+progress.show(len(tiers))+'\n')

    # primer searching in database by default
    if db:
        progress = Progressbar(len(intervals),'Querying database')
        for i, iv in enumerate(intervals):
            sys.stderr.write('\r'+progress.show(i))
            ivpairs[iv] = []
            primerpairs = db.query(iv)
            if primerpairs:
                for pair in primerpairs:
                    ivpairs[iv].append(pair)
                # remove excess primers (ordered by midpointdistance)
                ivpairs[iv] = ivpairs[iv][:config['report']['pairs']]
        sys.stderr.write('\r'+progress.show(len(intervals))+'\n')
        # print query count
        print >> sys.stderr, 'Found primers for {:d} out of {:d} intervals in database'.format(len([ iv for iv in intervals if ivpairs[iv]]), len(intervals))

    # designing
    if design:
        for tier in tiers:
            # get intervals which do not satisfy minimum amplicon number
            insufficentAmpliconIntervals = [ iv for iv in intervals if config['report']['pairs']>len(ivpairs[iv]) ]
            if not insufficentAmpliconIntervals:
                break  # end design process
            print >> sys.stderr, "Round #{} ({} intervals)".format(tier+1, len(insufficentAmpliconIntervals))
            # Primer3 design
            designedPairs = {}
            progress = Progressbar(len(insufficentAmpliconIntervals),'Designing primers')
            for i,iv in enumerate(insufficentAmpliconIntervals):
                sys.stderr.write('\r'+progress.show(i))
                try:
                    designIntervalOversize = max([ max(x) for x in config['design']['primer3'][tier]['PRIMER_PRODUCT_SIZE_RANGE'] ])
                except:
                    print >> sys.stderr, "WARNING: could not determine maximum amplicon size, default setting applied"
                    designIntervalOversize = 2000
                p3 = Primer3(config['design']['genome'], iv.locus(), designIntervalOversize)
                p3.design(iv.name, config['design']['primer3'][tier])
                if p3.pairs:
                    designedPairs[iv] = p3.pairs
                #else: print >> sys.stderr, '\n' +'\n'.join(p3.explain)
            sys.stderr.write('\r'+progress.show(len(insufficentAmpliconIntervals))+'\n')
            if designedPairs:
                ## import designed primer pairs (place on genome and get amplicons)
                with tempfile.NamedTemporaryFile(suffix='.fa',prefix="primers_",delete=False) as fh:
                    for k,v in designedPairs.items():
                        for pairnumber, pair in enumerate(v):
                            print >> fh, pair[0].fasta('_'.join([ k.name, str(pairnumber), 'rev' if k.strand < 0 else 'fwd' ]))
                            print >> fh, pair[1].fasta('_'.join([ k.name, str(pairnumber), 'fwd' if k.strand < 0 else 'rev' ]))
                pairs = importPrimerPairs(fh.name, config, primer3=True)
                os.unlink(fh.name)  # remove fasta file
                ## Remove non-specific and blacklisted primer pairs
                specificPrimerPairs = []
                blacklisted = 0
                for i, pair in enumerate(pairs):
                    if pair.uniqueid() in blacklist:
                        blacklisted += 1
                    elif all([pair[0].checkTarget(), pair[1].checkTarget()]):
                        specificPrimerPairs.append(pair)
                if blacklisted:
                    sys.stderr.write('Removed '+str(blacklisted)+' blacklisted primer pairs\n')
                if len(pairs) != len(specificPrimerPairs):
                    sys.stderr.write('Removed '+str(len(pairs)-len(specificPrimerPairs))+' non-specific primer pairs\n')
                pairs = specificPrimerPairs

                ## add SNPinfo (SNPcheck) for main target
                progress = Progressbar(len(pairs),'SNPcheck')
                for i, pair in enumerate(pairs):
                    sys.stderr.write('\r'+progress.show(i))
                    for p in pair:
                        p.snpCheckPrimer(config['snpcheck']['common'])
                sys.stderr.write('\r'+progress.show(len(pairs))+'\n')

                # assign designed primer pairs to intervals (remove ranks and tag)
                intervalindex = { iv.name: iv for iv in intervals }
                intervalprimers = { iv.name: set([ p.uniqueid() for p in ivpairs[iv] ]) for iv in intervals }
                failCount = 0
                for pair in pairs:
                    passed = 0
                    if pair.uniqueid() not in intervalprimers[pair.name]:
                        if pair.check(config['designlimits']):
                            # add default tag
                            for primer in pair:
                                primer.tag = config['design']['tag']
                            # assign to interval
                            ivpairs[intervalindex[pair.name]].append(pair)
                            intervalprimers[pair.name].add(pair.uniqueid())
                            # rename (for variant based naming which is too rich)
                            if rename:
                                ivpairs[intervalindex[pair.name]][-1].rename(rename)
                        else:
                            # add to blacklist if design limits fail
                            blacklist.append(pair.uniqueid())
                            failCount += 1
                if failCount:
                    print >> sys.stderr, 'INFO: {:2d} pairs violated design limits and were blacklisted ({} total)'.format(failCount,str(len(blacklist)))
                # print failed primer designs
                for k,v in intervalprimers.items():
                    if len(v)==0:
                        print >> sys.stderr, 'WARNING: Target {} failed on designlimits'.format(k)

    # save blacklist cache
    try:
        pickle.dump(list(set(blacklist)),open(config['blacklistcache'],'wb'))
    except:
        print >> sys.stderr, 'Could not write to blacklist cache, check permissions'
        print >> sys.stderr, os.getcwd(), config['blacklistcache']

    # print primer pair count and build database table
    failure = [ iv.name for iv,p in ivpairs.items() if config['report']['pairs']>len(p) ]
    print >> sys.stderr, 'got primers for {:d} out of {:d} targets'.format(len(ivpairs)-len(failure), len(ivpairs))
    print >> sys.stderr, '{:<20} {:9} {:<10}'.format('INTERVAL', 'AMPLICONS', 'STATUS')
    print >> sys.stderr, '-'*41
    for iv,p in sorted(ivpairs.items(),key=lambda x:x[0].name):
        print >> sys.stderr, '{:<20} {:9} {:<10}'.format(unquote(iv.name), len(p), "FAIL" if len(p)<config['report']['pairs'] else "OK")

    # select primer pairs
    primerTable = []  # primer table (text)
    primerVariants = defaultdict(list)  # primerpair -> intervalnames/variants dict
    missedIntervals = []  # list of missed intervals/variants
    print >> sys.stderr, '========'
    if compatible:
        # make pairs and exclude all pairings not in compatibility list score by geometric mean of 1based rank
        rankedPairs = []
        for l, left in enumerate(sorted(ivpairs[intervals[0]])):
            for r, rite in enumerate(sorted(ivpairs[intervals[1]])):
                if seqhash(left[0].seq, rite[1].seq) in compatible:
                    rankedPairs.append((((l+1)*(r+1)) ** (1.0/2), max(l, r), left, rite))
        if rankedPairs:
            rankedPairs.sort()  # sort ranked paired primerpairs
            primerTable.append([unquote(intervals[0].name)] + str(rankedPairs[0][2]).split('\t'))
            primerTable.append([unquote(intervals[1].name)] + str(rankedPairs[0][3]).split('\t'))
            primerVariants[rankedPairs[0][2]].append(intervals[0])
            primerVariants[rankedPairs[0][3]].append(intervals[1])
            # add gap interval
            gapInterval = Interval(intervals[0].chrom, intervals[0].chromStart, intervals[1].chromEnd)
            gapPrimerPair = PrimerPair([rankedPairs[0][2][0], rankedPairs[0][3][1]], name='gap_'+gapInterval.name )
            primerTable.append([unquote(gapInterval.name)] + str(gapPrimerPair).split('\t'))
        else:
            missedIntervals = intervals
    else:
        # select by best pairs independently (always print database primers)
        for iv in sorted(ivpairs.keys()):
            print "IV", unquote(iv.name)
            if not ivpairs[iv]:
                missedIntervals.append(iv)
            for i, p in enumerate(sorted(ivpairs[iv])):
                if i == config['report']['pairs']:
                    break  # only report number of primer pairs requested
                # log primer design (0 if from database)
                if p.designrank() >= 0:
                    p.log(config['logfile'])
                # save result (with interval names)
                primerVariants[p].append(iv)
                # save to primer table
                primerTable.append([unquote(iv.name)] + str(p).split('\t'))
    # update primer pairs with covered variants
    for pp, v in primerVariants.items():
        pp.variants = v
    return primerTable, primerVariants.keys(), missedIntervals

# ==============================================================================
# === convenience functions for webservice =====================================
# ==============================================================================

# query database / design primer for VCF,BED,GenePred or interval
def zippyPrimerQuery(config, targets, design=True, outfile=None, db=None, store=False, tiers=[0], gap=None):
    intervals = readTargets(targets, config['tiling'])  # get intervals from file or commandline
    if gap:  # gap PCR primers
        try:
            assert len(intervals)==1
            intervals += readTargets(gap, config['tiling'])  # get interval of second breakpoint
            assert len(set([i.chrom for i in intervals] ))
        except AssertionError:
            print >> sys.stderr, "ERROR: gap-PCR primers can only be designed for a single pair of breakpoint intervals on the same chromosome!"
        except:
            raise
    primerTable, resultList, missedIntervals = getPrimers(intervals,db,design,config,tiers,compatible=True if gap else False)
    ## print primerTable
    if outfile:
        with open(outfile,'w') as fh:
            print >> fh, '\n'.join([ '\t'.join(map(str,l)) for l in primerTable ])
    else:
        print >> sys.stdout, '\n'.join([ '\t'.join(map(str,l)) for l in primerTable ])
    ## print and store primer pairs
    # if db:
    if store and db and design:
        db.addPair(*resultList)  # store pairs in database (assume they are correctly designed as mispriming is ignored and capped at 1000)
        print >> sys.stderr, "Primer designs stored in database"
    return primerTable, resultList, missedIntervals

# batch query primer database and create confirmation worksheet
def zippyBatchQuery(config, targets, design=True, outfile=None, db=None, predesign=False, tiers=[0]):
    # read targets from first file and additional files
    if not isinstance(targets,list):
        targets = [ targets ]
    # load query files
    print >> sys.stderr, 'Reading batch file {}...'.format(targets[0])
    sampleVariants, genes, fullgenes = readBatch(targets[0], config['tiling'], database=db)
    for t in range(1,len(targets)): # read additional files
        print >> sys.stderr, 'Reading additional file {}...'.format(targets[t])
        sv, g, f = readBatch(targets[t], config['tiling'], database=db)
        # amend target regions
        for k,v in sv.items():
            if k in sampleVariants.keys():
                sampleVariants[k] += v
            else:
                sampleVariants[k] = v
        genes = list(set(genes) | set(g))
        fullgenes = list(set(fullgenes) | set(f))

    print >> sys.stderr, '\n'.join([ '{:<20} {:>2d}'.format(sample,len(variants)) \
        for sample,variants in sorted(sampleVariants.items(),key=lambda x: x[0]) ])

    # predesign
    if predesign and db and genes:
        designVariants = [ var for var in variants if not db.query(var) ]
        selectedgeneexons = list(set(genes) - set(fullgenes))
        print >> sys.stderr, "Designing exon primers for {} variants..".format(str(len(designVariants)))
        # get variants with no overlapping amplicon -> get variants which need new primer designs
        intervals = IntervalList([],source='GenePred')
        if designVariants:
            with open(config['design']['annotation']) as fh:
                for iv in GenePred(fh,getgenes=selectedgeneexons,**config['tiling']):  # get intervals from file or commandline
                    found = False
                    for dv in designVariants:
                        if not found and iv.overlap(dv):
                            intervals.append(iv)
                            found = True
                        if found: break
                    if found: break
        # add full genes
        if fullgenes:
            with open(config['design']['annotation']) as fh:
                intervals += GenePred(fh,getgenes=fullgenes,**config['tiling'])
        # predesign and store
        if intervals:
            primerTable, resultList, missedIntervals = getPrimers(intervals,db,predesign,config,tiers)
            if db:
                db.addPair(*resultList)  # store pairs in database (assume they are correctly designed as mispriming is ignored and capped at 1000)
        # reload query files ()
        print >> sys.stderr, 'Updating query table...'
        sampleVariants = readBatch(targets[0], config['tiling'], database=db)[0]
        for t in range(1,len(targets)): # read additional files
            sv = readBatch(targets[t], config['tiling'], database=db)[0]
            # amend target regions
            for k,v in sv.items():
                if k in sampleVariants.keys():
                    sampleVariants[k] += v
                else:
                    sampleVariants[k] = v

    # for each sample design
    primerTableConcat = []
    allMissedIntervals = {}
    missedIntervalNames = []
    tests = []  # tests to run
    for sample, intervals in sorted(sampleVariants.items(),key=lambda x: x[0]):
        print >> sys.stderr, "Getting primers for {} variants in sample {}".format(len(intervals),sample)
        # get/design primers
        #print >> sys.stderr, intervals
        primerTable, resultList, missedIntervals = getPrimers(intervals,db,design,config,tiers,rename=shortHumanReadable)
        if missedIntervals:
            allMissedIntervals[sample] = missedIntervals
            missedIntervalNames += [ i.name for i in missedIntervals ]
        # store result list
        primerTableConcat += [ [sample]+l for l in primerTable ]
        # store primers
        if db:
            db.addPair(*resultList)  # store pairs in database (assume they are correctly designed as mispriming is ignored and capped at 1000)
        # Build Tests
        for primerpair in resultList:
            tests.append(Test(primerpair,sample))
    ## print primerTable
    writtenFiles = []
    if not outfile:
        print >> sys.stdout, '\n'.join([ '\t'.join(l) for l in primerTableConcat ])
    else:
        worksheetName = '' if os.path.basename(targets[0]).startswith(os.path.basename(outfile)) else os.path.basename(outfile)
        # output data
        writtenFiles.append(outfile+'.txt')
        print >> sys.stderr, "Writing results to {}...".format(writtenFiles[-1])
        with open(writtenFiles[-1],'w') as fh:
            print >> fh, '\n'.join([ '\t'.join(l) for l in primerTableConcat ])
        # Primer Test Worksheet
        primerTests = [ t for t in tests if not any(t.primerpairobject.locations()) ]
        if primerTests:
            writtenFiles.append(outfile+'.primertest.pdf')
            print >> sys.stderr, "Writing Test Worksheet to {}...".format(writtenFiles[-1])
            ws = Worksheet(primerTests,name="Primer Test PCR")  # load worksheet
            ws.addControls(control='Normal')  # add positive controls
            ws.fillPlates(size=config['report']['platesize'],randomize=True,\
                includeSamples=False, includeControls=True)  # only include controls
            ws.createWorkSheet(writtenFiles[-1], primertest=True, worklist=worksheetName, **config['report'])
            # robot csv
            writtenFiles.append(outfile+'.primertest.csv')
            print >> sys.stderr, "Writing Test CSV to {}...".format(writtenFiles[-1])
            ws.robotCsv(writtenFiles[-1], sep=',')
            # order list
            writtenFiles.append(outfile+'.ordersheet.csv')
            print >> sys.stderr, "Writing primer order list to {}...".format(writtenFiles[-1])
            ws.orderCsv(writtenFiles[-1], config=config['ordersheet'])
        # Batch PCR worksheet
        writtenFiles.append(outfile+'.pdf')
        print >> sys.stderr, "Writing worksheet to {}...".format(writtenFiles[-1])
        ws = Worksheet(tests,name='Variant Confirmations')  # load worksheet
        ws.addControls()  # add controls
        ws.fillPlates(size=config['report']['platesize'],randomize=True)
        ws.createWorkSheet(writtenFiles[-1], worklist=worksheetName, **config['report'])
        # validate primer tube labels (checks for hash substring collisions)
        ws.tubeLabels()
        # robot csv
        writtenFiles.append(outfile+'.csv')
        print >> sys.stderr, "Writing robot CSV to {}...".format(writtenFiles[-1])
        ws.robotCsv(writtenFiles[-1], sep=',')
        # tube labels
        writtenFiles.append(outfile+'.tubelabels.txt')
        print >> sys.stderr, "Writing tube labels to {}...".format(writtenFiles[-1])
        ws.tubeLabels(writtenFiles[-1],tags=config['ordersheet']['sequencetags'])
        # write missed intervals
        missedIntervalNames = []
        if allMissedIntervals:
            writtenFiles.append(outfile+'.failed.txt')
            print >> sys.stderr, "Writing failed designs {}...".format(writtenFiles[-1])
            with open(writtenFiles[-1],'w') as fh:
                print >> fh, '\t'.join(['sample','variant'])
                for sample, missed in sorted(allMissedIntervals.items()):
                    print >> fh, '\n'.join([ '\t'.join([sample,unquote(i.name)]) for i in missed ])
    return writtenFiles, sorted(list(set(missedIntervalNames)))

# update storage location for primer
def updateLocation(primername, location, database, force=False):
    occupied = database.getLocation(location)
    if not occupied or force:
        if database.storePrimer(primername,location,force):
            print >> sys.stderr, '%s location sucessfully set to %s' % (primername, str(location))
            return ('success', location)
        else:
            print >> sys.stderr, 'WARNING: %s location update to %s failed' % (primername, str(location))
            return ('fail', location)
    else:
        print >> sys.stderr, 'Location already occupied by %s' % (' and '.join(occupied))
        return ('occupied', occupied)

# search primer pair by name substring matching
def searchByName(searchName, db):
    primersInDB = db.query(searchName)
    print >> sys.stderr, 'Found {} primer pairs with string "{}"'.format(len(primersInDB),searchName)
    return primersInDB

# update name of primer in database
def updatePrimerName(primerName, newName, db):
    nameUpdate = db.updateName(primerName, newName)
    if nameUpdate:
        print >> sys.stderr, 'Primer %s renamed %s' % (primerName, newName)
        return nameUpdate
    else:
        print >> sys.stderr, 'Primer renaming failed'
        return nameUpdate

# update name of primer pair in database
def updatePrimerPairName(pairName, newName, db):
    nameUpdate = db.updatePairName(pairName, newName)
    if nameUpdate:
        print >> sys.stderr, 'Pair %s renamed %s' % (pairName, newName)
        return nameUpdate
    else:
        print >> sys.stderr, 'Pair renaming failed'
        return nameUpdate

# blacklist primer pair in database
def blacklistPair(pairname, db):
    blacklisted = db.blacklist(pairname)
    print >> sys.stderr, '%s added to blacklist' % (blacklisted,)
    return blacklisted

# just delete primer pair in database (skip addition to blacklist)
def deletePair(pairname, db):
    deleted = db.blacklist(pairname,True)
    print >> sys.stderr, '%s has been deleted' % (deleted,)
    return deleted

def readprimerlocations(locationfile):
    header = []
    updateList = []
    with open(locationfile) as csvfile:
        readfile = csv.reader(csvfile, delimiter=',')
        for line in readfile:
            if not header:
                header = line
            elif len(line)<3:
                pass  # empty line (CR/LF)
            else:
                try:
                    row = dict(zip(header,line))
                    updateList.append([row['PrimerName'], Location(row['Box'].strip('Bbox'), row['Well'])])
                except:
                    raise
                    raise Exception('InputFormatError')
    return updateList

# ==============================================================================
# === CLI ======================================================================
# ==============================================================================
def main():
    from zippylib import ascii_encode_dict
    from zippylib import banner

    print >> sys.stderr, banner(__version__)

    parser = ArgumentParser(prog="zippy.py", description= 'Zippy - Primer design and database')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__+'('+__status__+')',\
        help="Displays version")

    #   configuration files
    global_group = parser.add_argument_group('Global options')
    global_group.add_argument("-c", dest="config", default='zippy.json',metavar="JSON_FILE", \
        help="configuration file [zippy.json]")
    global_group.add_argument("--tiers", dest="tiers", default='0,1,2', \
        help="Allowed design tiers (0,1,...,n)")

    # run modes
    subparsers = parser.add_subparsers(help='Help for subcommand')

    ## add primers
    parser_add = subparsers.add_parser('add', help='Add previously designed primers to database')
    parser_add.add_argument("primers", default=None, metavar="FASTA/TAB", \
        help="Primers or locations to add to database")
    parser_add.set_defaults(which='add')

    ## retrieve
    parser_retrieve = subparsers.add_parser('get', help='Get/design primers')
    parser_retrieve.add_argument("targets", default=None, metavar="VCF/BED/Interval/GenePred", \
        help="File with intervals of interest or CHR:START-END (mandatory for gap-PCR)")
    parser_retrieve.add_argument("--design", dest="design", default=False, action="store_true", \
        help="Design primers if not in database")
    parser_retrieve.add_argument("--gap", dest="gap", default=None, metavar="CHR:START-END", \
        help="Second break point for gap-PCR")
    parser_retrieve.add_argument("--nostore", dest="store", default=True, action='store_false', \
        help="Do not store result in database")
    parser_retrieve.add_argument("--outfile", dest="outfile", default='', type=str, \
        help="Output file name")
    parser_retrieve.set_defaults(which='get')

    ## query database for primers by name
    parser_query = subparsers.add_parser('query', help='Query database for primers with specified sub-string in name')
    parser_query.add_argument("subString", default=None, metavar="Sub-string within name", \
        help="String found within primer name")
    parser_query.set_defaults(which='query')

    ## batch
    parser_batch = subparsers.add_parser('batch', help='Batch design primers for sample list')
    parser_batch.add_argument("targets", default=None, metavar="FILE1,FILE2,...", \
        help="SNPpy result table(s) ")
    parser_batch.add_argument("--predesign", dest="predesign", default=False, action="store_true", \
        help="Design primers for all genes in batch")
    parser_batch.add_argument("--nodesign", dest="design", default=True, action="store_false", \
        help="Skip primer design if not in database")
    parser_batch.add_argument("--tiers", dest="tiers", default='0,1,2', \
        help="Allowed design tiers (0,1,...,n)")
    parser_batch.add_argument("--outfile", dest="outfile", default='', type=str, \
        help="Create worksheet PDF, order and robot CSV")
    parser_batch.set_defaults(which='batch')

    ## update
    parser_update = subparsers.add_parser('update', help='Update status and location of primers')
    parser_update.add_argument('-l', dest="location", nargs=3, \
        help="Update storage location of primer pair (primerid vessel well)")
    parser_update.add_argument('-t', dest="locationtable", \
        help="Batch update storage locations from TSV (primerid vessel well)")
    parser_update.add_argument("--force", dest="force", default=False, action='store_true', \
        help="Force Location update (resets existing)")
    parser_update.add_argument('-b', dest="blacklist", type=str, \
        help="Blacklist primer")
    parser_update.set_defaults(which='update')

    ## dump specific datasets from database
    parser_dump = subparsers.add_parser('dump', help='Data dump')
    parser_dump.add_argument("--amplicons", dest="amplicons", default='', type=str, \
        help="Retrieve amplicons of given size (eg. 10-1000)")
    parser_dump.add_argument("--ordersheet", dest="ordersheet", default=False, action="store_true", \
        help="IDT order sheet (primer pairs with no status marker)")
    parser_dump.add_argument("--locations", dest="locations", default=False, action="store_true", \
        help="Primer locations")
    parser_dump.add_argument("--redundancies", dest="redundancies", default=False, action="store_true", \
        help="Primers with same sequence and tag")
    parser_dump.add_argument("--table", dest="table", default=False, action="store_true", \
        help="Primer pair table with locations")
    parser_dump.add_argument("--outfile", dest="outfile", default='', type=str, \
        help="Output file name")
    parser_dump.set_defaults(which='dump')

    options = parser.parse_args()

    # read config and open database
    with open(options.config) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
    here = config['primerbed'] if 'primerbed' in config.keys() and config['primerbed'] else None
    db = PrimerDB(config['database'],dump=here)

    if options.which=='add':  # read primers and add to database
        # import primer pairs
        if options.primers.split('.')[-1].startswith('fa'):
            pairs = importPrimerPairs(options.primers, config, primer3=False)  # import and locate primer pairs
            print >> sys.stderr, "Storing Primers..."
            db.addPair(*pairs)  # store pairs in database (assume they are correctly designed as mispriming is ignored and capped at 1000)
            sys.stderr.write('Added {} primer pairs to database\n'.format(len(pairs)))
        # store locations if table
        if not options.primers.split('.')[-1].startswith('fa'):  # assume table format
            locations = importPrimerLocations(options.primers)
            print >> sys.stderr, "Setting Primer locations..."
            db.addLocations(*locations.items())
            sys.stderr.write('Added {} locations for imported primers\n'.format(len(locations)))
    elif options.which=='dump':  # data dump fucntions (`for bulk downloads`)
        if options.amplicons:
            try:
                l = options.amplicons.split('-')
                assert len(l)==2
                amplen = map(int,l)
            except (AssertionError, ValueError):
                raise ConfigError('must give amplicon size to retrieve')
            except:
                raise
            else:
                # get amplicons amplen
                data,colnames = db.dump('amplicons', size=amplen)
        elif options.ordersheet:
            data,colnames = db.dump('ordersheet', **config['ordersheet'])
        elif options.locations:
            data,colnames = db.dump('locations')
        elif options.table:
            data,colnames = db.dump('table')
        elif options.redundancies:
            data,colnames = db.getRedundantPrimers()
        else:
            print >> sys.stderr, "What to dump stranger?"
            sys.exit(1)
        # format data output
        if options.outfile:
            dump = Data(data,colnames)
            dump.writefile(options.outfile)  # sets format by file extension
        else:
            print '\t'.join(colnames)
            for row in data:
                print '\t'.join(map(str,row))
    elif options.which=='update':  #update location primer pairs are stored
        if options.location:
            primer, vessel, well = options.location
            updateLocation(primer, Location(vessel, well), db, options.force)
        if options.locationtable:
            updateList = readprimerlocations(options.locationtable)
            for item in updateList:  # [ Primer, Location ]
                updateLocation(item[0], item[1], db, options.force)
        if options.blacklist:
            print >> sys.stderr, 'BLACKLISTED PAIRS: {}'.format(','.join(db.blacklist(options.blacklist)))
            print >> sys.stderr, 'REMOVED ORPHANS:   {}'.format(','.join(db.removeOrphans()))
    elif options.which=='get':  # get primers for targets (BED/VCF or interval)
        zippyPrimerQuery(config, options.targets, options.design, options.outfile, \
            db, options.store, map(int,options.tiers.split(',')), options.gap)
    elif options.which=='batch':
        zippyBatchQuery(config, options.targets.split(','), options.design, options.outfile, \
            db, options.predesign, map(int,options.tiers.split(',')))
    elif options.which=='query':
        searchByName(options.subString, db)

if __name__=="__main__":
    main()
