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
from primulalib.vcf import VCF
from primulalib.primer import MultiFasta, Primer3
from primulalib.database import PrimerDB
from argparse import ArgumentParser


# functionalities
'''
import primers from sequence
    blat
    insert into db

design primers from Vcf
    parse Vcf
    check database
    design primer pairs
'''

parser = ArgumentParser(prog="primula.py", \
    description= 'Primula - Primer design and database')
parser.add_argument('--version', action='version', version='%(prog)s '+__version__+'('+__status__+')',\
    help="Displays version")

#   configuration files
config_group = parser.add_argument_group('Configuration options')
config_group.add_argument("-c", dest="config", default='primula.json',metavar="JSON_FILE", \
    help="configuration file [primula.json]")

#   pipeline options
pipeline_group = parser.add_argument_group('Run options')
pipeline_group.add_argument("-p", "--primers", dest="primers", default=None, metavar="FASTA/TAB", \
    help="Primer FASTA to add to database (automatically finds targets)")
pipeline_group.add_argument("-v", "--vcf", dest="vcf", default=None, metavar="VCF", \
    help="Variant file with variants to design")

# Run options
run_group = parser.add_argument_group('Run options')
run_group.add_argument("--database", dest="database", default=False, action="store_true", \
    help="Query database")
run_group.add_argument("--design", dest="design", default=False, action="store_true", \
    help="Design primers if not in database")

options = parser.parse_args()

## read configuration
with open(options.config) as conf:
    config = json.load(conf)

# open database handler
db = PrimerDB(config['database'])

# read primers and add to database
if options.primers:
    primerfile = MultiFasta(options.primers)
    primers = primerfile.createPrimers(config['import']['bowtieindex'])  # places in genome
    for primer in primers:
        print >> sys.stderr, "Adding primer", primer
        primer.calcProperties()  # calc Tm and GC
        db.addPrimer(primer)

    print repr(db)

elif options.vcf:
    with open(options.vcf) as fh:
        vcf = VCF(fh)
    for variant in vcf:
        print variant
        if options.database:
            primerpair = db.query(variant, config)  # get primer pair based on Tm and distance (and variant size)
            if primerpair:
                pass
                #output
        if options.design:
            p3 = Primer3(config['primer3']['genome'])
            #design primer
            primerpair = p3.design('test',variant.locus(),config['primer3']['settings'])
            print primerpair
            sys.exit('debugStop')
            #STORE?

# change stock?
# elif options.stock:
#     pass
#     # change stock
