# Zippy
Primer database and design tool

## Description
This program integrates a simple SQLite primer database and design tool based on the Primer3 library.
It allows the automatic generation of primer pairs based on a VCF, BED or SNPpy result table.

### Primer design
Zippy's primary functionality is the automatic design of working sequencing primer pairs. It uses Primer3 as a backend and verifies the uniqueness of resulting sequening amplicons and identifies overlap with potentially interfering common variation in the genome. Multiple design parameter sets can be specified and Zippy will autmatically move on to the next, less stringent parameter sets if the design process does not yield any good primer pairs (deep searching).

### Storage Management
Zippy also manages storage of primer pairs in boxes (vessel) and wells. It also provides functionality for barcode tracking of primer dilutions ready for sequencing.

### Batch processing
Zippy encourages a workflow of variant confirmations from Variant scoring with SNPpy. It directly reads output from such and generates multiple outputs:

1. Overview of sequencing primers selected for variant confirmation (TXT)
2. Worksheet to assist and track the confirmation workflow (PDF)
3. Program instruction for Hamilton robot to prepare plates for Sanger sequencing confirmations (CVS)
4. Barcoded tube labels for Primer dilutions (ZPL)
5. Worksheet and Program for test of new primers (optional,PDF/CSV)
6. A list of variants with no suitable sequencing primers (optional,TXT)

For information on the Hamilton Robot program please contact the authors.

## Install
The installation procedure installs an Apache2 webserver and installs all required python modules in a *virtualenv* in `/usr/local/zippy`. Genomic data is stored

### Virtual machine setup
fire up the virtual machine and connect with
> `vagrant up && vagrant ssh`

then follow the local installation instructions below.

### Webservice Install
The current installation routine will download and build the human *GRCh_37* genome index and download the common variantion data from *dbsnp142*. If you desire to use your own resources or an alternative reference genome simply put everything into the `./resource` directory and it will be imported to `/var/local/zippy/resources` during the installation process.
Make sure to modify the configuration file `zippy.json` accordingly.

The webservice will be exposed to the VM host at address `55.55.55.5`.

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

Currently the design process is executed synchronously. This can potentially lead to timeouts during the design process if many target regions are requested. In this case please run *zippy* from the command line interface. This will be changed in a future version.

### Command line interface
Before running zippy from the CLI, make sure to activate the virtual environment first
> `source /usr/local/zippy/venv/bin/activate`

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
- Primer3 base design
- genome mispriming check (bowtie)
- SNPcheck (crossvalidate primer location with indexed VCF)
- amplicon size validation on import
- batch processing of SNPpy results

### v1.1
- worksheet and robot CSV file generation

### v1.2
- Webinterface/GUI
- Tube Label Generation
- Fixed Plate filling
- GenePred input for design

### v2.0
- New database schema to allow for primer collections
- Primer tag tracking
- Primer table import
- Storage by Primer and storage validation
- Allows for same primer in multiple pairs
- added primer redundancy query
- Primer design parameter sets (deep digging for primers)
- Importing of primer pairs with multiple amplicons
- Blacklist cache for design
- Improved webinterface
- Apache Webserver in VM
- New setup routines (zippyprimer on PyPi)
- Easier VM provisioning

### FUTURE
- Support for primer collections (multiplexing)
- Asynchronous design process on webinterface
- Web GUI extensions
- Storage map (suggest new locations?)
- Import from files (fasta,list)

## Reference Genome
Primers are designed from the unmasked 1kg reference genome while ignoring simple repeat regions.
Alternatively, use a Repeatmasked reference genome. Even better if common SNPs are masked as well (eg. >1%).
BEDTOOLS offers `maskfasta` for this purpose. Masking data in BED format can be obtained from UCSC (http://genome.ucsc.edu/cgi-bin/hgTables?command=start).
