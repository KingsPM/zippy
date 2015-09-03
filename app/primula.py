#!/usr/bin/env python

__doc__=="""
################################################################
# Primula - Primer database and automated design               #
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

import sys
import json
from primulalib.files import VCF, BED, Interval
from primulalib.primer import MultiFasta, Primer3, Primer
from primulalib.database import PrimerDB
from argparse import ArgumentParser
from collections import defaultdict

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
    left_suffix, rite_suffix = ['l','L','5'],['r','R','3']
    pairs = []
    for i in range(len(primers)):
        primerI = primers[i].name.split('_')
        for j in range(i,len(primers)):
            primerJ = primers[j].name.split('_')
            if '_'.join(primerI[:-1]) == '_'.join(primerJ[:-1]):
                if primerI[-1][0] in left_suffix and primerJ[-1][0] in rite_suffix:
                    pairs.append([primers[i], primers[j]])
                elif primerJ[-1][0] in left_suffix and primerI[-1][0] in rite_suffix:
                    pairs.append([primers[j], primers[i]])
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
    parser = ArgumentParser(prog="primula.py", description= 'Primula - Primer design and database')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__+'('+__status__+')',\
        help="Displays version")

    #   configuration files
    config_group = parser.add_argument_group('Configuration options')
    config_group.add_argument("-c", dest="config", default='primula.json',metavar="JSON_FILE", \
        help="configuration file [primula.json]")

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
        db.addPair(*pairs)  # store pairs in database (assume they are correctly designed as mispriming is ignored)
        print repr(db)  # show database content (debugging only)
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
        for iv in intervals:
            print iv
            if options.design and iv not in ivpairs.keys():  # not in database

                p3 = Primer3(config['primer3']['genome'],iv.locus(),1000)  # genome and target
                p3.design('test',config['primer3']['settings'])
                print p3.pairs[0][0].meta
                p3.show()  # show placed primers
                sys.exit()

                print '\n'.join([ str(i)+':'+str(v) for i,v in enumerate(p3.pairs)])
                continue

                for pair in p3.pairs:
                    print pair[0].name, pair[1].name
                ivpairs[iv] = p3.pairs

        # check new designs for mispriming and create valid primer pairs
        # for k,v in ivpairs:
        #     print pair[0].name, pair[1].name


    # change stock?
    # elif options.stock:
    #     pass
    #     # change stock
