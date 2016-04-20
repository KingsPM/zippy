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
__version__ = "2.0.0"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import os
import re
import sys
import json
import tempfile
import hashlib
from zippylib.files import VCF, BED, GenePred, Interval, Data, readTargets, readBatch
from zippylib.primer import Genome, MultiFasta, Primer3, Primer, PrimerPair, Location, parsePrimerName
from zippylib.reports import Test
from zippylib.database import PrimerDB
from zippylib import ConfigError, Progressbar, banner
from zippylib.reports import Worksheet
from argparse import ArgumentParser
from collections import defaultdict, Counter
import cPickle as pickle

blacklistCacheFile = '.blacklist.cache'

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
                    l['vessel'] and l['well']:
                    # store location
                    loc = Location(l['vessel'],l['well'])
                    if l['primername'] in primerlocs.keys():
                        primerlocs[l['primername']].merge(loc)
                    else:
                        primerlocs[l['primername']] = loc
    return primerlocs

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
    primers = primerfile.createPrimers(config['design']['bowtieindex'],delete=False,tags=primertags)  # places in genome
    # pair primers (by name or by primerset)
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
                    pairs[setname][1] = p
                else:
                    if pairs[setname][0] is None:
                        pairs[setname][0] = p
                    else:
                        assert pairs[setname][1] is None
                        pairs[setname][1] = p
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
        acceptedPairs = []
        print >> sys.stderr, 'Identifying correct amplicons for unplaced primer pairs...'
        for p in pairs.values():
            if not p[0].targetposition or not p[1].targetposition:
                amplicons = p.amplicons(config['import']['ampliconsize'],autoreverse=True)
                if amplicons:
                    shortest = sorted(amplicons,key=lambda x: len(x[2]))[0]  # sort amplicons by size
                    if len(amplicons)>1:
                        print >> sys.stderr, 'WARNING: mutiple amplicons for {}. Assuming shortest ({}bp) is correct.'.format(p.name,str(len(shortest[2])))
                    p[0].targetposition = shortest[0]  # m
                    p[1].targetposition = shortest[1]  # n
                    acceptedPairs.append(p)
                elif not primer3:  # try to find amplicon by sequence matching
                    refGenome = Genome(config['design']['genome'])
                    if p[0].loci and len(p[0].loci) < len(p[1].loci):
                        query, mapped = 1,0
                    elif p[1].loci and len(p[1].loci) < len(p[0].loci):
                        query, mapped = 0,1
                    else:
                        print >> sys.stderr, 'WARNING: Primer set {} is too ambiguous to be imported ({},{})'.format(p.name,len(p[0].loci),len(p[1].loci))
                        continue
                    # find by quick sequence matching
                    for l in p[mapped].loci:
                        loc = refGenome.primerMatch(l,p[query].seq,config['import']['ampliconsize'])
                        if loc:
                            p[query].loci += loc
                    # store new amplicon
                    amplicons = p.amplicons(config['import']['ampliconsize'],autoreverse=True)
                    if amplicons:
                        p[0].targetposition = amplicons[0][0]  # m
                        p[1].targetposition = amplicons[0][1]  # n
                        acceptedPairs.append(p)
                    else:
                        print >> sys.stderr, '\n'.join([ ">"+str(l) for l in p[0].loci])
                        print >> sys.stderr, '\n'.join([ "<"+str(l) for l in p[1].loci])
                        print >> sys.stderr, 'WARNING: Primer set {} has no valid amplicons ({},{})'.format(p.name,len(p[0].loci),len(p[1].loci))
                else:
                    print >> sys.stderr, 'WARNING: Primer set {} does not produce a well-sized, unique amplicon ({},{})'.format(p.name,len(p[0].loci),len(p[1].loci))
            else:
                acceptedPairs.append(p)
        validPairs = acceptedPairs
    return validPairs

'''get primers from intervals'''
def getPrimers(intervals, db, design, config, deep=True):
    ivpairs = defaultdict(list)  # found/designed primer pairs (from database or design)
    blacklist = db.blacklist() if db else []
    try:
        blacklist += pickle.load(open(blacklistCacheFile,'rb'))
    except:
        print >> sys.stderr, 'Could not read blacklist cache, check permissions'
        print >> sys.stderr, os.getcwd(), blacklistCacheFile

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
        maxTier = len(config['design']['primer3']) if deep else 1  # only search first tier unless deep
        for tier in range(maxTier):
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
        pickle.dump(list(set(blacklist)),open(blacklistCacheFile,'wb'))
    except:
        print >> sys.stderr, 'Could not write to blacklist cache, check permissions'
        print >> sys.stderr, os.getcwd(), blacklistCacheFile

    # print primer pair count and build database table
    failure = [ iv.name for iv,p in ivpairs.items() if config['report']['pairs']>len(p) ]
    print >> sys.stderr, 'got primers for {:d} out of {:d} targets'.format(len(ivpairs)-len(failure), len(ivpairs))
    print >> sys.stderr, '{:<20} {:9} {:<10}'.format('INTERVAL', 'AMPLICONS', 'STATUS')
    print >> sys.stderr, '-'*41
    for iv,p in sorted(ivpairs.items(),key=lambda x:x[0].name):
        print >> sys.stderr, '{:<20} {:9} {:<10}'.format(iv.name, len(p), "FAIL" if len(p)<config['report']['pairs'] else "OK")

    ## get best primer pairs
    ##### PRIORITISE AND ALWAYS PRINT DATABASE PRIMERS
    print >> sys.stderr, '========'
    primerTable, resultList, missedIntervals = [], [], []
    for iv in sorted(ivpairs.keys()):
        if not ivpairs[iv]:
            missedIntervals.append(iv)
        for i, p in enumerate(sorted(ivpairs[iv])):
            if i == config['report']['pairs']:
                break  # only report number of primer pairs requested
            resultList.append(p)
            if p.designrank() >= 0:
                p.log(config['logfile'])
            primerTable.append([iv.name] + str(p).split('\t'))
    return primerTable, resultList, missedIntervals

# ==============================================================================
# === convenience functions ====================================================
# ==============================================================================

def zippyPrimerQuery(config, targets, design=True, outfile=None, db=None, store=False, deep=True):
    intervals = readTargets(targets, config['tiling'])  # get intervals from file or commandline
    # get/design primer pairs
    primerTable, resultList, missedIntervals = getPrimers(intervals,db,design,config,deep)
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

def zippyBatchQuery(config, targets, design=True, outfile=None, db=None, predesign=False, deep=True):
    print >> sys.stderr, 'Reading batch file {}...'.format(targets)
    sampleVariants, genes = readBatch(targets, config['tiling'])
    print >> sys.stderr, '\n'.join([ '{:<20} {:>2d}'.format(sample,len(variants)) \
        for sample,variants in sorted(sampleVariants.items(),key=lambda x: x[0]) ])
    # predesign
    if predesign and genes and db:
        print >> sys.stderr, "Designing primers for {} genes..".format(str(len(genes)))
        # parse gene intervals from refGene
        with open(config['design']['annotation']) as fh:
            intervals = GenePred(fh,getgenes=genes,**config['tiling'])  # get intervals from file or commandline
        # predesign and store
        primerTable, resultList, missedIntervals = getPrimers(intervals,db,predesign,config,deep)
        if db:
            db.addPair(*resultList)  # store pairs in database (assume they are correctly designed as mispriming is ignored and capped at 1000)
    # for each sample design
    primerTableConcat = []
    allMissedIntervals = {}
    tests = []  # tests to run
    for sample, intervals in sorted(sampleVariants.items(),key=lambda x: x[0]):
        print >> sys.stderr, "Getting primers for {} variants in sample {}".format(len(intervals),sample)
        # get/design primers
        primerTable, resultList, missedIntervals = getPrimers(intervals,db,design,config,deep)
        if missedIntervals:
            allMissedIntervals[sample] = missedIntervals
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
            ws.createWorkSheet(writtenFiles[-1], primertest=True, **config['report'])
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
        ws = Worksheet(tests,name='Validation batch PCR')  # load worksheet
        ws.addControls()  # add controls
        ws.fillPlates(size=config['report']['platesize'],randomize=True)
        ws.createWorkSheet(writtenFiles[-1],**config['report'])
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
                    print >> fh, '\n'.join([ '\t'.join([sample,i.name]) for i in missed ])
                    missedIntervalNames += [ i.name for i in missed ]
    return writtenFiles, sorted(list(set(missedIntervalNames)))

def updateLocation(primername, location, database):
    occupied = database.getLocation(location)
    if not occupied:
        if database.storePrimer(primername,location):
            return '%s location sucessfully set to %s' % (primername, str(location))
        else:
            return 'WARNING: %s location update to %s failed' % (primername, str(location))
    else:
        return 'Location already occupied by %s' % (' and '.join(occupied))

def searchByName(searchName, db):
    primersInDB = db.queryName(searchName)
    print >> sys.stderr, 'Found {} primer pairs with string "{}"'.format(len(primersInDB),searchName)
    return primersInDB

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

    # run modes
    subparsers = parser.add_subparsers(help='Help for subcommand')

    ## add primers
    parser_add = subparsers.add_parser('add', help='Add previously designed primers to database')
    parser_add.add_argument("primers", default=None, metavar="FASTA/TAB", \
        help="Primer FASTA to add to database (automatically finds targets)")
    parser_add.set_defaults(which='add')

    ## retrieve
    parser_retrieve = subparsers.add_parser('get', help='Get/design primers')
    parser_retrieve.add_argument("targets", default=None, metavar="VCF/BED/Interval/GenePred", \
        help="File with intervals of interest or chr:start-end")
    parser_retrieve.add_argument("--design", dest="design", default=False, action="store_true", \
        help="Design primers if not in database")
    parser_retrieve.add_argument("--nodeep", dest="deep", default=True, action='store_false', \
        help="Skip deep search for primers")
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
    parser_batch.add_argument("targets", default=None, metavar="SNPpy result table", \
        help="SNPpy result table")
    parser_batch.add_argument("--predesign", dest="predesign", default=False, action="store_true", \
        help="Design primers for all genes in batch")
    parser_batch.add_argument("--nodesign", dest="design", default=True, action="store_false", \
        help="Skip primer design if not in database")
    parser_batch.add_argument("--nodeep", dest="deep", default=True, action='store_false', \
        help="Skip deep search for primers")
    parser_batch.add_argument("--outfile", dest="outfile", default='', type=str, \
        help="Create worksheet PDF, order and robot CSV")
    parser_batch.set_defaults(which='batch')

    ## update
    parser_update = subparsers.add_parser('update', help='Update status and location of primers')
    parser_update.add_argument('-l', dest="location", nargs=3, \
        help="Update storage location of primer pair (pairid vessel well)")
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
    parser_dump.add_argument("--outfile", dest="outfile", default='', type=str, \
        help="Output file name")
    parser_dump.set_defaults(which='dump')

    options = parser.parse_args()

    # read config and open database
    with open(options.config) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
    db = PrimerDB(config['database'])

    if options.which=='add':  # read primers and add to database
        # import primer pairs
        pairs = importPrimerPairs(options.primers, config, primer3=False)  # import and locate primer pairs
        db.addPair(*pairs)  # store pairs in database (assume they are correctly designed as mispriming is ignored and capped at 1000)
        sys.stderr.write('Added {} primer pairs to database\n'.format(len(pairs)))
        # store locations if table
        if not options.primers.split('.')[-1].startswith('fa'):  # assume table format
            locations = importPrimerLocations(options.primers)
            db.addLocations(*locations.items())
            sys.stderr.write('Added {:2d} locations for imported primers\n'.format(len(locations)))
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
            print >> sys.stderr, updateLocation(primer, Location(vessel, well), db)
        if options.blacklist:
            print >> sys.stderr, 'BLACKLISTED PAIRS: {}'.format(','.join(db.blacklist(options.blacklist)))
            print >> sys.stderr, 'REMOVED ORPHANS:   {}'.format(','.join(db.removeOrphans()))
    elif options.which=='get':  # get primers for targets (BED/VCF or interval)
        zippyPrimerQuery(config, options.targets, options.design, options.outfile, db, options.store, options.deep)
    elif options.which=='batch':
        zippyBatchQuery(config, options.targets, options.design, options.outfile, db, options.predesign, options.deep)
    elif options.which=='query':
        searchByName(options.subString, db)

if __name__=="__main__":
    main()
