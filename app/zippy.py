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
__version__ = "1.0"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import os
import re
import sys
import json
import tempfile
from zippylib.files import VCF, BED, Interval, Data, readTargets, readBatch
from zippylib.primer import MultiFasta, Primer3, Primer, PrimerPair
from zippylib.database import PrimerDB
from zippylib import ConfigError, Progressbar, banner

from argparse import ArgumentParser
from collections import defaultdict, Counter

'''reads fasta fastafile, searches genome for targets, decides if valid pairs'''
def importPrimerPairs(fastafile,primer3=True):
    # read primers (fasta file)
    primerfile = MultiFasta(fastafile)
    print >> sys.stderr, "Placing primers on genome..."
    # Align primers to genome and add Tm/GC
    primers = primerfile.createPrimers(config['targeting']['bowtieindex'])  # places in genome
    # add Tm/GC
    for primer in primers:
        primer.calcProperties()
    print >> sys.stderr, "Calculated Tm/GC for {} primers".format(len(primers))
    # pair primers
    left_suffix, rite_suffix = ['F','f','L','l','5','left'],['R','r','3','right']
    pairs = {}
    for p in primers:
        # get name and suffix
        primername, primersuffix = '_'.join(p.name.split('_')[:-1]), p.name.split('_')[-1]
        try:
            pairs[primername]
        except KeyError:
            pairs[primername] = PrimerPair([None,None])
        except:
            raise
        if primersuffix in left_suffix:
            pairs[primername][0] = p
        elif primersuffix in rite_suffix:
            pairs[primername][1] = p
        else:
            raise Exception('PrimerNameError')

    # check if any unpaired primers
    validPairs = [ p for p in pairs.values() if all(p) ]
    assert len(validPairs) == len(pairs.values())

    if primer3:  # prune ranks and read target
        for p in validPairs:
            p.pruneRanks()
    else:  # guess target if not set
        acceptedPairs = []
        print >> sys.stderr, 'Identifying correct amplicons for unplaced primer pairs...'
        for p in validPairs:
            if not p[0].targetposition or not p[1].targetposition:
                if len(p.amplicons(config['import']['ampliconsize']))==1 or len(p.reverse().amplicons(config['import']['ampliconsize']))==1:
                    amplicons = p.amplicons(config['import']['ampliconsize'])
                    p[0].targetposition = amplicons[0][0]  # m
                    p[1].targetposition = amplicons[0][1]  # n
                    acceptedPairs.append(p)
                else:
                    if p.reversed:
                        p.reverse()
                    print >> sys.stderr, 'WARNING: Primer {} does not produce a single, well-sized amplicon ({},{})'.format(p.name(),len(p[0].loci),len(p[1].loci))
                    #print >> sys.stderr, '\tFWD', ','.join([ str(l) for l in p[0].loci ])
                    #print >> sys.stderr, '\tREV', ','.join([ str(l) for l in p[1].loci ])
                    continue
            else:
                acceptedPairs.append(p)
        validPairs = acceptedPairs
    return validPairs

'''get primers from intervals'''
def getPrimers(intervals, options):
    ivpairs = defaultdict(list)  # found/designed primer pairs (from database or design)
    blacklist = db.blacklist()
    # primer searching in database by default
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
    if options.debug:
        for iv,pp in ivpairs.items():
            if not pp:
                print >> sys.stderr, 'no primer for', iv
    print >> sys.stderr, 'Found primers for {:d} out of {:d} intervals in database'.format(len([ iv for iv in intervals if ivpairs[iv]]), len(intervals))

    # show blacklist
    if options.debug:
        bl = db.blacklist()
        print >> sys.stderr, '\n++BLACKLIST+++++++++++++++++++++++++++'
        print >> sys.stderr, '\n'.join(bl) if bl else 'empty'
        print >> sys.stderr, '++++++++++++++++++++++++++++++++++++++\n'

    # designing
    if options.design:
        designedPairs = {}
        progress = Progressbar(len(intervals),'Designing primers')
        for i,iv in enumerate(intervals):
            sys.stderr.write('\r'+progress.show(i))
            if iv not in ivpairs.keys() or config['report']['pairs']>len(ivpairs[iv]):  # not in database or not enough primer pairs for interval
                p3 = Primer3(config['primer3']['genome'], iv.locus(), 300)  # genome and target region (plusminus)
                p3.design(iv.name, config['primer3']['settings'])
                if p3.pairs:
                    designedPairs[iv] = p3.pairs
                else:
                    print >> sys.stderr, '\n' +'\n'.join(p3.explain)
        sys.stderr.write('\r'+progress.show(len(intervals))+'\n')

        if designedPairs:
            ## place primer pairs
            with tempfile.NamedTemporaryFile(suffix='.fa',prefix="primers_",delete=False) as fh:
                for k,v in designedPairs.items():
                    for pairnumber, pair in enumerate(v):
                        print >> fh, pair[0].fasta('_'.join([ k.name, str(pairnumber), "left" ]))
                        print >> fh, pair[1].fasta('_'.join([ k.name, str(pairnumber), "right" ]))
            pairs = importPrimerPairs(fh.name,primer3=True)
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
            checkedPairs = []  # checked primer pairs (with correct target)
            for i, pair in enumerate(pairs):
                sys.stderr.write('\r'+progress.show(i))
                for p in pair:
                    p.snpCheckPrimer(config['snpcheck']['common'])
            sys.stderr.write('\r'+progress.show(len(pairs))+'\n')

            # assign designed primer pairs to intervals (and remove ranks)
            intervalindex = { iv.name: iv for iv in intervals }
            intervalprimers = { iv.name: set([ p.uniqueid() for p in ivpairs[iv] ]) for iv in intervals }
            for pair in pairs:
                passed = 0
                if pair.uniqueid() not in intervalprimers[pair.name()]:
                    if pair.check(config['designlimits']):
                        ivpairs[intervalindex[pair.name()]].append(pair)
                        intervalprimers[pair.name()].add(pair.uniqueid())
            # print failed primer designs
            for k,v in intervalprimers.items():
                if len(v)==0:
                    print >> sys.stderr, 'Primer {} failed on designlimits'.format(k)

    # print primer pair count and build database table
    failure = [ iv.name for iv,p in ivpairs.items() if config['report']['pairs']>len(p) ]
    print >> sys.stderr, 'got primers for {:d} out of {:d} targets'.format(len(ivpairs)-len(failure), len(ivpairs))
    print >> sys.stderr, '{:<20} {:9} {:<10}'.format('INTERVAL', 'AMPLICONS', 'STATUS')
    print >> sys.stderr, '-'*41
    for iv,p in sorted(ivpairs.items(),key=lambda x:x[0].name):
        print >> sys.stderr, '{:<20} {:9} {:<10}'.format(iv.name, len(p), "!!" if len(p)<config['report']['pairs'] else "")

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
            primerTable.append([iv.name] + repr(p).split())

    return primerTable, resultList, missedIntervals


''' main script '''
if __name__=="__main__":

    # print banner
    print >> sys.stderr, banner(__version__)

    parser = ArgumentParser(prog="zippy.py", description= 'Zippy - Primer design and database')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__+'('+__status__+')',\
        help="Displays version")

    #   configuration files
    global_group = parser.add_argument_group('Global options')
    global_group.add_argument("-c", dest="config", default='zippy.json',metavar="JSON_FILE", \
        help="configuration file [zippy.json]")
    global_group.add_argument("--debug", dest="debug", default=False, action="store_true", \
        help="Debugging (show database dump at end)")

    # run modes
    subparsers = parser.add_subparsers(help='help for subcommand')

    ## add primers
    parser_add = subparsers.add_parser('add', help='Add previously designed primers to database')
    parser_add.add_argument("primers", default=None, metavar="FASTA/TAB", \
        help="Primer FASTA to add to database (automatically finds targets)")
    parser_add.set_defaults(which='add')

    ## retrieve
    parser_retrieve = subparsers.add_parser('get', help='Get/design primers')
    parser_retrieve.add_argument("targets", default=None, metavar="VCF/BED/Interval", \
        help="File with intervals of interest or chr:start-end")
    parser_retrieve.add_argument("--design", dest="design", default=False, action="store_true", \
        help="Design primers if not in database")
    parser_retrieve.add_argument("--nostore", dest="nostore", default=False, action='store_true', \
        help="Do not store result in database")
    parser_retrieve.add_argument("--outfile", dest="outfile", default='', type=str, \
        help="Output file name")
    parser_retrieve.set_defaults(which='get')

    ## batch
    parser_batch = subparsers.add_parser('batch', help='Batch design primers for sample list')
    parser_batch.add_argument("targets", default=None, metavar="SNPpy result table", \
        help="SNPpy result table")
    parser_batch.add_argument("--design", dest="design", default=True, action="store_true", \
        help="Design primers if not in database [TRUE]")
    parser_batch.add_argument("--outfile", dest="outfile", default='', type=str, \
        help="Output base name (order, worksheet, unavailable)")
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
    parser_dump.add_argument("--amplicons", dest="amplicons", default='10-1000', type=str, \
        help="Retrieve amplicons of given size [10-1000]")
    parser_dump.add_argument("--ordersheet", dest="ordersheet", default=False, action="store_true", \
        help="IDT order sheet (primer pairs with no status marker)")
    parser_dump.add_argument("--outfile", dest="outfile", default='', type=str, \
        help="Output file name")
    parser_dump.set_defaults(which='dump')

    options = parser.parse_args()

    ## read configuration (convert unicode to ascii string)
    def ascii_encode_dict(data):
        ascii_encode = lambda x: x.encode('ascii') if type(x) is unicode else x
        return dict(map(ascii_encode, pair) for pair in data.items())
    with open(options.config) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)

    # open database handler
    db = PrimerDB(config['database'])

    if options.which=='add':  # read primers and add to database
        pairs = importPrimerPairs(options.primers, primer3=False)  # import and locate primer pairs
        db.addPair(*pairs)  # store pairs in database (assume they are correctly designed as mispriming is ignored and capped at 1000)
        sys.stderr.write('Added {} primer pairs to database\n'.format(len(pairs)))
        if options.debug: print repr(db)  # show database content (debugging only)
    elif options.which=='dump':  # data dump fucntions (`for bulk downloads`)
        # dump amplicons fo given size to stdout
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
            if options.ordersheet:
                data,colnames = db.dump('ordersheet',size=amplen, **config['ordersheet'])
            else:
                data,colnames = db.dump('amplicons',size=amplen)
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
            pairid = options.location[0]
            vessel = options.location[1]
            well = options.location[2]
            if not db.storePrimer(pairid,vessel,well):
                print >> sys.stderr, 'Location already occupied' # Try and include statement of primer pair stored at location
            else:
                print >> sys.stderr, 'Primer pair location updated'
        if options.blacklist:
            db.blacklist(options.blacklist)
    elif options.which=='get':  # get primers for targets (BED/VCF or interval)
        intervals = readTargets(options.targets, config['tiling'])  # get intervals from file or commandline

        # get/design primer pairs
        primerTable, resultList, missedIntervals = getPrimers(intervals,options)

        ## print primerTable
        if options.outfile:
            with open(options.outfile,'w') as fh:
                print >> fh, '\n'.join([ '\t'.join(l) for l in primerTable ])
        else:
            print >> sys.stdout, '\n'.join([ '\t'.join(l) for l in primerTable ])

        ## print and store primer pairs
        if not options.nostore:
            db.addPair(*resultList)  # store pairs in database (assume they are correctly designed as mispriming is ignored and capped at 1000)

        ## print database content
        if options.debug:
            print >> sys.stderr, '\n++DATABASE++++++++++++++++++++++++++++'
            print >> sys.stderr, repr(db)
            print >> sys.stderr, '++++++++++++++++++++++++++++++++++++++\n'

    elif options.which=='batch':
        sampleVariants = readBatch(options.targets, config['tiling'])
        print >> sys.stderr, '\n'.join([ '{:<20} {:>2d}'.format(sample,len(variants)) \
            for sample,variants in sorted(sampleVariants.items(),key=lambda x: x[0]) ])
        # for each sample design
        primerTableConcat = []
        for sample, intervals in sorted(sampleVariants.items(),key=lambda x: x[0]):
            print >> sys.stderr, "Getting primers for {} variants in sample {}".format(len(intervals),sample)
            # get/design primers
            primerTable, resultList, missedIntervals = getPrimers(intervals,options)

            # store result list
            primerTableConcat += [ [sample]+l for l in primerTable ]

            # store primers
            db.addPair(*resultList)  # store pairs in database (assume they are correctly designed as mispriming is ignored and capped at 1000)

        ## print primerTable
        print >> sys.stdout, '\n'.join([ '\t'.join(l) for l in primerTableConcat ])

        ## reports
        # add controls

        # primer ordering

        # fill plate sheet

        # print ordersheet
