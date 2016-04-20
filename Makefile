# installation makefile

ZIPPYPATH=/usr/local/zippy
ZIPPYVAR=/var/local/zippy
ZIPPYWWW=/var/www

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

zippy-install:
	mkdir -p $(ZIPPYPATH)
	cd $(ZIPPYPATH) && virtualenv venv
	$(ZIPPYPATH)/venv/bin/pip install cython
	$(ZIPPYPATH)/venv/bin/pip install pysam
	$(ZIPPYPATH)/venv/bin/pip install primer3-py
	$(ZIPPYPATH)/venv/bin/pip install -r zippy/bpackage-requirements.txt


apache: apache-user apache-install apache-restart

apache-install:
	apt-get install -y apache2 apache2.2-common apache2-mpm-prefork apache2-utils libexpat1 ssl-cert
	apt-get install -y libapache2-mod-wsgi

apache-zippy:
	rsync -a zippy.wsgi $(ZIPPYWWW)
	rsync -a resources $(ZIPPYVAR)
	cp apache2_zippy /etc/apache2/sites-available/zippy
	sudo a2ensite zippy

apache-user:
	useradd -M flask
	usermod -s /bin/false flask
	usermod -L flask
	adduser flask www-data
	chown -R flask:www-data /var/local/zippy
	chown -R flask:www-data /var/www/zippy

apache-restart: apache-zippy
	/etc/init.d/apache2 restart

bowtie:
	wget -c http://netix.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip && \
	unzip bowtie2-2.2.6-linux-x86_64.zip && \
	cd bowtie2-2.2.6 && \
	mv bowtie2* /usr/local/bin
	rm -rf bowtie2-2.2.6

genome: genome-download genome-index

genome-download:
	mkdir -p $(ZIPPYVAR) && cd $(ZIPPYVAR) && \
	wget ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz && \
	wget ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai && \
	gunzip human_g1k_v37.fasta.gz

genome-index:
	cd $(ZIPPYVAR) && \
	bowtie2-build human_g1k_v37.fasta human_g1k_v37.bowtie

annotation: variation-download refgene-download

variation-download:
	echo "Please provide VCF file for SNP check (eg. dbsnp_138.b37.G5A.vcf.gz)"
	echo "AND COPY TO" $(ZIPPYVAR)

refgene-download:
	mkdir -p $(ZIPPYVAR) && cd $(ZIPPYVAR) && \
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -D hg19 -P 3306 \
	 -e 'SELECT DISTINCT r.bin,CONCAT(r.name,".",i.version),c.ensembl,r.strand,
	r.txStart,r.txEnd,r.cdsStart,r.cdsEnd,r.exonCount,r.exonStarts,r.exonEnds,
	r.score,r.name2,r.cdsStartStat,r.cdsEndStat,r.exonFrames
	FROM refGene as r, gbCdnaInfo as i, ucscToEnsembl as c
	WHERE r.name=i.acc AND c.ucsc = r.chrom
	ORDER BY r.bin;' > refGene
