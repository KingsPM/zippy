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
__version__ = "0.9"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Development"

import os
import re
import sys
import json
import tempfile
from zippylib.files import VCF, BED, Interval, Data, readTargets
from zippylib.primer import MultiFasta, Primer3, Primer, PrimerPair
from zippylib.database import PrimerDB
from zippylib import ConfigError, Progressbar, banner

from argparse import ArgumentParser
from collections import defaultdict, Counter

'''
reads fasta fastafile
searches genome for targets
decides if valid pairs
'''
def importPrimerPairs(fastafile):
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
    # remove unpaired
    validPairs = [ p for p in pairs.values() if all(p) ]
    print >> sys.stderr, "Got %s paired primers" % str(len(validPairs))
    return validPairs  # return valid pairs

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
        help="Debugging")
    global_group.add_argument("--quiet", dest="quiet", default=False, action="store_true", \
        help="Minimise screen output (dont print results)")
    global_group.add_argument("--outfile", dest="outfile", default='', type=str, \
        help="Output file name (bed,interval,fasta)")

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
    parser_retrieve.set_defaults(which='get')

    ## dump specific datasets from database
    parser_dump = subparsers.add_parser('dump', help='Data dump')
    parser_dump.add_argument("--amplicons", dest="amplicons", default='', type=str, \
        help="Retrieve possible amplicons of given size (eg. 100-800)")
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
        pairs = importPrimerPairs(options.primers)  # import (and locate) primer pairs
        db.addPair(*pairs)  # store pairs in database (assume they are correctly designed as mispriming is ignored and capped at 1000)
        sys.stderr.write('Added {} primer pairs to database\n'.format(len(pairs)))
        if options.debug: print repr(db)  # show database content (debugging only)
    elif options.which=='dump':  # data dump fucntions (`for bulk downloads`)
        if options.amplicons:
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
                data,colnames = db.dump('amplicons',size=amplen)

            # format data output
            if options.outfile:
                dump = Data(data,colnames)
                dump.writefile(options.outfile)  # sets format by file extension
            else:
                print '\t'.join(colnames)
                for row in data:
                    print '\t'.join(map(str,row))
    elif options.which=='get':  # get primers for targets (BED/VCF or interval)
        intervals = readTargets(options.targets, config['tiling'])  # get intervals from file or commandline
        ivpairs = defaultdict(list)  # found/designed primer pairs (from database or design)
        blacklist = set()
        # primer searching in database by default
        progress = Progressbar(len(intervals),'Querying database')
        for i, iv in enumerate(intervals):
            sys.stderr.write('\r'+progress.show(i))
            primerpairs = db.query(iv, config['tiling']['flank'])
            for pair in primerpairs:
                if pair.status is None or status != 0:
                    ivpairs[iv].append(pair)
                else:
                    blacklist.add(hash(pair))
        sys.stderr.write('\r'+progress.show(len(intervals))+'\n')
        print >> sys.stderr, 'Found primers for {:d} out of {:d} intervals in database'.format(len([ iv for iv in intervals if ivpairs[iv]]), len(intervals))

        # designing
        if options.design:
            designedPairs = {}
            progress = Progressbar(len(intervals),'Designing primers')
            for i,iv in enumerate(intervals):
                sys.stderr.write('\r'+progress.show(i))
                if options.debug:
                    print iv
                if iv not in ivpairs.keys() or config['report']['pairs']>len(ivpairs[iv]):  # not in database or not enough primer pairs for interval
                    p3 = Primer3(config['primer3']['genome'],iv.locus(),300)  # genome and target
                    p3.design(iv.name,config['primer3']['settings'])
                    if options.debug:
                        print >> sys.stderr
                        p3.show()  # show placed primers
                        print '\n'.join([ str(i)+':'+str(v) for i,v in enumerate(p3.pairs)])
                        for pair in p3.pairs:
                            print pair[0].name, pair[1].name
                    designedPairs[iv] = p3.pairs
            sys.stderr.write('\r'+progress.show(len(intervals))+'\n')

            if designedPairs:
                ## place primer pairs
                with tempfile.NamedTemporaryFile(suffix='.fa',prefix="primers_",delete=False) as fh:
                    for k,v in designedPairs.items():
                        for pairnumber, pair in enumerate(v):
                            print >> fh, pair[0].fasta('_'.join([ k.name, str(pairnumber), "left" ]))
                            print >> fh, pair[1].fasta('_'.join([ k.name, str(pairnumber), "right" ]))
                pairs = importPrimerPairs(fh.name)
                os.unlink(fh.name)  # remove fasta file

                ## Remove non-specific and blacklisted primer pairs
                specificPrimers = []
                blacklisted = 0
                for i, pair in enumerate(pairs):
                    if hash(pair) in blacklist:
                        blacklisted += 1
                    elif all([pair[0].checkTarget(), pair[1].checkTarget()]):
                        specificPrimers.append(pair)
                if blacklisted:
                    sys.stderr.write('Removed '+str(blacklisted)+' blacklisted primer pairs\n')
                if len(pairs) != len(specificPrimers):
                    sys.stderr.write('Removed '+str(len(pairs)-len(specificPrimers))+' non-specific primer pairs\n')
                pairs = specificPrimers

                ## add SNPinfo (SNPcheck) for main target
                progress = Progressbar(len(pairs),'SNPcheck')
                checkedPairs = []  # checked primer pairs (with correct target)
                for i, pair in enumerate(pairs):
                    sys.stderr.write('\r'+progress.show(i))
                    for p in pair:
                        p.snpCheckPrimer(config['snpcheck']['common'])
                sys.stderr.write('\r'+progress.show(len(pairs))+'\n')

                # assign designed primer pairs to intervals
                intervalindex = { i.name: i for i in intervals }
                for pair in pairs:
                    intervalName = '_'.join(pair[0].name.split('_')[:-2])
                    if intervalName not in ivpairs.keys():
                        ivpairs[intervalindex[intervalName]].append(pair)

            ## add ivpairs with no primers
            for iv in set(intervals).difference(set(ivpairs.keys())):
                ivpairs[iv] = []

        # print primer pair count and build database table
        failure = len([ iv.name for iv,p in ivpairs.items() if config['report']['pairs']>len(p) ])
        print >> sys.stderr, 'got primers for {:d} out of {:d} targets'.format(len(ivpairs)-failure, len(ivpairs))
        if options.quiet:
            print >> sys.stderr, "FAILED INTERVALS ({}): {}".format(len(failure), ','.join(failure))
        else:
            print >> sys.stderr, '{:<16} {:>3} {:<10}'.format('INTERVAL', 'AMP', 'STATUS')
            print >> sys.stderr, '-'*31
            for iv,p in sorted(ivpairs.items(),key=lambda x:x[0].name):
                print >> sys.stderr, '{:<16} {:>3} {:<10}'.format(iv.name, len(p), "!!" if len(p)<config['report']['pairs'] else "")

        ## get best primer pairs
        resultList = []
        for iv in ivpairs.keys():
            for i, p in enumerate(sorted(ivpairs[iv])):
                if i == config['report']['pairs']: break  # only report number of primer pairs requested
                #if False in sortvalues(p): continue  ##DATABASE DOES NOT RETURN ATTRIBUTES YET
                resultList.append(p)

        ## store primer pairs
        db.addPair(*resultList)  # store pairs in database (assume they are correctly designed as mispriming is ignored and capped at 1000)

        # dump database (debugging)
        # print >> sys.stderr, '++++++++++++++++++'
        # print >> sys.stderr, "DB DUMP:"
        # print >> sys.stderr, repr(db)
        # print >> sys.stderr, '++++++++++++++++++'

        # WRITE RESULT PRIMERS
        if options.outfile:
            fh = open(options.outfile,'w') if options.outfile != '-' else sys.stdout
            for pair in resultList:
                for p in pair:
                    p.strand = '-' if p.targetposition.reverse else '+'
                    print >> fh, ">",p.name+"|"+p.targetposition.chrom+":"+str(p.targetposition.offset)+\
                    "-"+str(int(p.targetposition.offset)+int(p.targetposition.length))+"\n"+\
                    str(p.seq)
            try:
                fh.close()
            except:
                pass

        # WRITE RESULT PRIMERS
        if not options.quiet:
            print >> sys.stderr, '==== RESULTS ===='
            print >> sys.stderr, '\n'.join([ repr(r) for r in resultList])

    # change stock?
    # elif options.stock:
    #     pass
    #     # change stock
