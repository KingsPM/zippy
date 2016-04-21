# installation makefile

ZIPPYPATH=/usr/local/zippy
ZIPPYVAR=/var/local/zippy
ZIPPYWWW=/var/www/zippy

WWWUSER=flask
WWWGROUP=www-data

all: install genome annotation

install: essential bowtie zippy-install apache

essential:
	locale-gen en_GB.UTF-8
	apt-get update
	apt-get install -y sqlite3 unzip git htop
	apt-get install -y python-pip python2.7-dev ncurses-dev python-virtualenv
	apt-get install -y libxslt-dev libxml2-dev libffi-dev
	apt-get install -y redis-server
	apt-get install -y build-essential libjpeg-dev libfreetype6-dev python-dev python-imaging libcurl3-dev
	apt-get install -y mysql-client

bowtie:
	wget -c http://netix.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip && \
	unzip bowtie2-2.2.6-linux-x86_64.zip && \
	cd bowtie2-2.2.6 && mv bowtie2* /usr/local/bin
	rm -rf bowtie2-2.2.6 bowtie2-2.2.6-linux-x86_64.zip

zippy-install:
	mkdir -p $(ZIPPYPATH)
	cd $(ZIPPYPATH) && virtualenv venv
	$(ZIPPYPATH)/venv/bin/pip install cython
	$(ZIPPYPATH)/venv/bin/pip install pysam
	$(ZIPPYPATH)/venv/bin/pip install primer3-py==0.5.0
	$(ZIPPYPATH)/venv/bin/pip install -r package-requirements.txt

apache: apache-install apache-user apache-setup apache-restart

apache-install:
	apt-get install -y apache2 apache2.2-common apache2-mpm-prefork apache2-utils libexpat1 ssl-cert
	apt-get install -y libapache2-mod-wsgi

apache-user:
	useradd -M $(WWWUSER)
	usermod -s /bin/false $(WWWUSER)
	usermod -L $(WWWUSER)
	adduser $(WWWUSER) $(WWWGROUP)

apache-setup: import-resources
	# make WWW directories
	mkdir -p $(ZIPPYWWW)
	rsync -a zippy.wsgi $(ZIPPYWWW)
	chown -R $(WWWUSER):$(WWWGROUP) $(ZIPPYWWW)
	# apache WSGI config
	cp apache2_zippy /etc/apache2/sites-available/zippy
	a2ensite zippy
	# disable default site
	a2dissite default

apache-restart:
	# restart apache
	/etc/init.d/apache2 restart

import-resources:
	# Copy resource files
	mkdir -p $(ZIPPYVAR)
	rsync -avPp resources $(ZIPPYVAR)
	chown -R $(WWWUSER):$(WWWGROUP) $(ZIPPYVAR)

genome: genome-download genome-index

genome-download:
	mkdir -p $(ZIPPYVAR)/resources && cd $(ZIPPYVAR)/resources && \
	wget ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz && \
	wget ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai && \
	gunzip human_g1k_v37.fasta.gz

genome-index:
	cd $(ZIPPYVAR)/resources && \
	bowtie2-build human_g1k_v37.fasta human_g1k_v37.bowtie

annotation: variation-download refgene-download

variation-download:
	mkdir -p $(ZIPPYVAR)/resources && cd $(ZIPPYVAR)/resources && \
	wget ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/common_all_20160408.vcf.gz && \
	wget ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/common_all_20160408.vcf.gz.tbi

refgene-download:
	mkdir -p $(ZIPPYVAR)/resources && cd $(ZIPPYVAR)/resources && \
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -D hg19 -P 3306 \
	 -e 'SELECT DISTINCT r.bin,CONCAT(r.name,".",i.version),c.ensembl,r.strand,
	r.txStart,r.txEnd,r.cdsStart,r.cdsEnd,r.exonCount,r.exonStarts,r.exonEnds,
	r.score,r.name2,r.cdsStartStat,r.cdsEndStat,r.exonFrames
	FROM refGene as r, gbCdnaInfo as i, ucscToEnsembl as c
	WHERE r.name=i.acc AND c.ucsc = r.chrom
	ORDER BY r.bin;' > refGene
