# Zippy
Primer database and design tool

## Description
This program integrates a simple SQLite primer database and design tool based on the Primer3 library.
It allows the automatic generation of primer pairs based on a VCF, BED or SNPpy result table.

## Install
### Virtual machine setup
fire up the virtual machine and connect with
> `vagrant up && vagrant ssh`

then follow the local installation instructions below.

### Webservice Install
The current installation routine will download and build the human *GRCh_37* genome index and download the common variantion data from *dbsnp142*. If you desire to use your own resources or an alternative reference genome simply put everything into the `./resource` directory and it will be imported to `/var/local/zippy/resources` during the installation process.
Make sure to modify the configuration file `zippy.json` accordingly.

You can install zippy with all needes resource with
> `sudo make all`

To install zippy and the flask webservice on apache2/wsgi_mod run
> `sudo make install`

Download b37 genomes and index with
> `sudo make genome`

Download variantion data and reference annotation with
> `sudo make annotation`

## Usage

### Webservice
The application runs on Apache Webserver (WSGI).
The standard install exposes the service on port 5000.
In the VM this port is mapped to 80 on the host machine.

### CLI
To design primers and query existing
> `zippy.py get <VCF/BED> --design`

add primers to database
> `zippy.py add <FASTA>`

retrieve (and design) primer for a single location
> `zippy.py get chr1:202030-202100 --design`

design/retrieve primers for a SNPpy result (batch mode)
> `zippy.py batch <SNPpy>`

Blacklist primer
> `zippy.py update -b <PRIMERNAME>`

Set/update primer storage location
> `zippy.py update -l <PRIMERNAME> <VESSEL> <WELL>`


## Release Notes
### v1.0
Primer3 base design
genome mispriming check (bowtie)
SNPcheck (crossvalidate primer location with indexed VCF)
amplicon size validation on import
batch processing of SNPpy results

### v1.1
worksheet and robot CSV file generation

### v1.2
Webinterface/GUI
Tube Label Generation
Fixed Plate filling
GenePred input for design

### v2.0
New database schema to allow for primer collections
Primer tag tracking
Primer table import
Storage by Primer and storage validation
Allows for same primer in multiple pairs
added primer redundancy query
Primer design parameter sets (deep digging for primers)
Importing of primer pairs with multiple amplicons
Blacklist cache for design
Improved webinterface
Apache Webserver in VM
New setup routines (zippyprimer on PyPi)
Easier VM provisioning

### FUTURE
Support for primer collections (multiplexing)
Asynchronous design process (celery)
Web GUI extensions
Storage map (suggest new locations?)
Import from files (fasta,list)

## Reference Genome
Primers are designed from the unmasked 1kg reference genome while ignoring simple repeat regions.
Alternatively, use a Repeatmasked reference genome. Even better if common SNPs are masked as well (eg. >1%).
BEDTOOLS offers `maskfasta` for this purpose. Masking data in BED format can be obtained from UCSC (http://genome.ucsc.edu/cgi-bin/hgTables?command=start).
