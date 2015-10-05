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
__version__ = "0.1"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Development"

import os
import sys
import json
import tempfile
from zippylib.files import VCF, BED, Interval, Data
from zippylib.primer import MultiFasta, Primer3, Primer
from zippylib.database import PrimerDB
from zippylib import ConfigError

from argparse import ArgumentParser
from collections import defaultdict, Counter

'''
reads fasta fastafile
searches genome for targets
decides if valid pairs
'''
def importPrimerPairs(fastafile):
    primerfile = MultiFasta(fastafile)
    primers = primerfile.createPrimers(config['targeting']['bowtieindex'])  # places in genome
    for primer in primers:
        primer.calcProperties()  # calc Tm and GC
    # find pairs
    left_suffix, rite_suffix = ['F','f','L','l','5'],['R','r','3']
    pairs = []
    for i in range(len(primers)):
        primerI = [ primers[i].name[:-1], primers[i].name[-1:] ]
        for j in range(i,len(primers)):
            if i==j: continue  # skip same primer
            primerJ = [ primers[j].name[:-1], primers[j].name[-1:] ]
            if primerI[0] == primerJ[0]:
                if primerI[1] in left_suffix and primerJ[1] in rite_suffix:
                    pairs.append([primers[i], primers[j]])
                elif primerJ[1] in left_suffix and primerI[1] in rite_suffix:
                    pairs.append([primers[j], primers[i]])
                else:
                    raise Exception('saasd')
    # return valid pairs
    return pairs

''' read target intervals from VCF, BED or directly'''
def readTargets(targets):
    with open(targets) as fh:
        if targets.endswith('vcf'):
            intervals = VCF(fh,**config['tiling'])
        elif targets.endswith('bed'):
            intervals = BED(fh,**config['tiling'])
        elif re.match('\w+:\d+-\d+',targets):
            m = re.match('(\w+):(\d+)-(\d+)',targets)
            intervals = [ Interval(*m.groups()) ]
        else:
            raise Exception('UnkownFile')
    return intervals


'''
import primers from sequence
    blat
    insert into db

design primers from Vcf
    parse Vcf
    check database
    design primer pairs
'''
if __name__=="__main__":
    parser = ArgumentParser(prog="zippy.py", description= 'Zippy - Primer design and database')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__+'('+__status__+')',\
        help="Displays version")

    #   configuration files
    config_group = parser.add_argument_group('Configuration options')
    config_group.add_argument("-c", dest="config", default='zippy.json',metavar="JSON_FILE", \
        help="configuration file [zippy.json]")
    config_group.add_argument("--debug", dest="debug", default=False, action="store_true", \
        help="Debugging")

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
    parser_retrieve.add_argument("--database", dest="database", default=False, action="store_true", \
        help="Query database")
    parser_retrieve.add_argument("--deep", dest="deep", default=False, action="store_true", \
        help="Allow new primer combinations")
    parser_retrieve.add_argument("--design", dest="design", default=False, action="store_true", \
        help="Design primers if not in database")
    parser_retrieve.set_defaults(which='get')

    ## dump specific datasets from database
    parser_dump = subparsers.add_parser('dump', help='Data dump')
    parser_dump.add_argument("--amplicons", dest="amplicons", default='', type=str, \
        help="Retrieve possible amplicons of given size (eg. 100-800)")
    parser_dump.add_argument("--outfile", dest="outfile", default='', type=str, \
        help="Output file name (bed,interval,fasta)")
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
        intervals = readTargets(options.targets)  # get intervals from file or commandline
        ivpairs = {}  # found/designed primer pairs (from database or design)
        # primer searching
        for iv in intervals:
            if options.database:  # check if inteval covered by primer pair
                raise NotImplementedError
                ivpairs[iv] = db.query(iv, config['primer3']['settings'])
                if ivpairs[iv]:
                    print "Found %d pairs for iv %s" % (len(ivpairs[iv]), iv)
                elif options.deep:  ## check if a new combination of primers would work
                    raise NotImplementedError
        # designing
        print >> sys.stderr, 'Designing primers'
        for i,iv in enumerate(intervals):
            print >> sys.stderr, '\r'+str(i)+'/'+str(len(intervals)),
            if options.debug:
                print iv
            if options.design and iv not in ivpairs.keys():  # not in database
                p3 = Primer3(config['primer3']['genome'],iv.locus(),300)  # genome and target
                p3.design(iv.name,config['primer3']['settings'])
                if options.debug:
                    p3.show()  # show placed primers
                if options.debug:
                    print '\n'.join([ str(i)+':'+str(v) for i,v in enumerate(p3.pairs)])
                    for pair in p3.pairs:
                        print pair[0].name, pair[1].name
                ivpairs[iv] = p3.pairs

        # check new designs for mispriming and create valid primer pairs
        ## write to fasta (process batch at once
        print >> sys.stderr, "\rChecking genome wide mispriming"
        fh = tempfile.NamedTemporaryFile(suffix='.fa',prefix="primers_",delete=False)
        for k,v in ivpairs.items():
            print k
            for pairnumber, pair in enumerate(v):
                print >> fh, pair[0].fasta('_'.join([ k.name, str(pairnumber), "left" ]))
                print >> fh, pair[1].fasta('_'.join([ k.name, str(pairnumber), "right" ]))
        fh.close()
        ## create primers with mispriming added
        pairs = importPrimerPairs(fh.name)
        ## remove fasta file
        os.unlink(fh.name)
        ## add SNPinfo (SNPcheck)
        counter = Counter()
        for p in pairs:
            counter['pairs'] += 1
            if len(p[0].loci)==1 and len(p[1].loci)==1:
                print p[0].loci[0]
                print p[1].loci[0]
                # is unique
                counter['unique'] += 1
                p[0].snpCheck(config['snpcheck']['common'])
                p[1].snpCheck(config['snpcheck']['common'])
                # no snp in the last third of the sequence (3 prime)
                if len(p[0].snp)==0 and len(p[1].snp==0):
                    counter['unique'] += 1
                    # no snp at ALL
                    print p[0],p[1]
        ## print primer



    # change stock?
    # elif options.stock:
    #     pass
    #     # change stock
