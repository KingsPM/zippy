# installation makefile

ZIPPYPATH=/usr/local/zippy
ZIPPYVAR=/var/local/zippy
ZIPPYWWW=/var/www/zippy

WWWUSER=flask
WWWGROUP=www-data

# production install
release: install resources webservice

# development installs (with mounted volume)
all: install resources

install: essential bowtie zippy-install

# requirements
essential:
	apt-get install -y wget
	apt-get install -y sqlite3 unzip git htop
	apt-get install -y python-pip python2.7-dev ncurses-dev python-virtualenv
	apt-get install -y libxslt-dev libxml2-dev libffi-dev
	apt-get install -y redis-server
	apt-get install -y build-essential libjpeg-dev libfreetype6-dev python-dev python-imaging libcurl3-dev
	apt-get install -y mysql-client
	apt-get install -y postgres-server-dev-9.1
	# add apache user
	useradd -M $(WWWUSER)
	usermod -s /bin/false $(WWWUSER)
	usermod -L $(WWWUSER)
	adduser $(WWWUSER) $(WWWGROUP)
	# install apache/wsgi
	apt-get install -y apache2 apache2.2-common apache2-mpm-prefork apache2-utils libexpat1 ssl-cert
	apt-get install -y libapache2-mod-wsgi
	# disable default site
	a2dissite 000-default

bowtie:
	wget -c https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.3.1/bowtie2-2.3.3.1-linux-x86_64.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbowtie-bio%2Ffiles%2Fbowtie2%2F2.3.3.1&ts=1507890163&use_mirror=kent && \
	unzip bowtie2-2.3.3.1-linux-x86_64.zip && \
	cd bowtie2-2.3.3.1-linux-x86_64 && mv bowtie2* /usr/local/bin
	rm -rf bowtie2-2.3.3.1*

# zippy setup (will move to distutils in future release)
zippy-install:
	# virtualenv
	mkdir -p $(ZIPPYPATH)
	cd $(ZIPPYPATH) && virtualenv venv
	$(ZIPPYPATH)/venv/bin/pip install Cython==0.24
	$(ZIPPYPATH)/venv/bin/pip install -r package-requirements.txt
	# create empty database
	mkdir -p $(ZIPPYVAR)
	touch $(ZIPPYVAR)/zippy.sqlite
	touch $(ZIPPYVAR)/zippy.log
	touch $(ZIPPYVAR)/.blacklist.cache
	mkdir -p $(ZIPPYVAR)/uploads
	mkdir -p $(ZIPPYVAR)/results
	chmod -R 777 $(ZIPPYVAR)

# gunicorn/nginx webserver
unicorn:
	apt-get install nginx
	# start gunicorn with
	# gunicorn --bind 0.0.0.0:8000 wsgi:app

# webservice install (production)
webservice:
	rsync -a --exclude-from=.gitignore . $(ZIPPYPATH)
	# make WWW directories
	mkdir -p $(ZIPPYWWW)
	cp install/zippy.wsgi $(ZIPPYWWW)/zippy.wsgi
	chown -R $(WWWUSER):$(WWWGROUP) $(ZIPPYWWW)
	# apache WSGI config
	cp install/zippy.hostconfig /etc/apache2/sites-available/zippy.conf
	# enable site and restart
	a2ensite zippy
	/etc/init.d/apache2 restart

# same for development environment (not maintained)
webservice-dev:
	# make WWW directories
	mkdir -p $(ZIPPYWWW)
	cp install/zippy_dev.wsgi $(ZIPPYWWW)/zippy.wsgi
	chown -R $(WWWUSER):$(WWWGROUP) $(ZIPPYWWW)
	# apache WSGI config
	cp install/zippy_dev.hostconfig /etc/apache2/sites-available/zippy.conf
	# enable site and restart
	a2ensite zippy
	/etc/init.d/apache2 restart

#### genome resources
import-resources:
	# Copy resource files
	mkdir -p $(ZIPPYVAR)/resources
	rsync -avPp resources $(ZIPPYVAR)
	chown -R $(WWWUSER):$(WWWGROUP) $(ZIPPYVAR)

resources: genome annotation

genome: genome-download genome-index

genome-download:
	mkdir -p $(ZIPPYVAR)/resources
	cd $(ZIPPYVAR)/resources && \
	wget -qO- ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz | \
	gzip -dcq | cat > human_g1k_v37.fasta && rm -f human_g1k_v37.fasta.gz && \
	wget -c ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai

genome-index:
	cd $(ZIPPYVAR)/resources && \
	bowtie2-build human_g1k_v37.fasta human_g1k_v37.bowtie

annotation: variation-download refgene-download

variation-download:
	mkdir -p $(ZIPPYVAR)/resources && cd $(ZIPPYVAR)/resources && \
	wget -c ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/00-common_all.vcf.gz && \
	wget -c ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/00-common_all.vcf.gz.tbi

refgene-download:
	mkdir -p $(ZIPPYVAR)/resources && cd $(ZIPPYVAR)/resources && \
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -D hg19 -P 3306 \
	 -e "SELECT DISTINCT r.bin,CONCAT(r.name,'.',i.version),c.ensembl,r.strand, r.txStart,r.txEnd,r.cdsStart,r.cdsEnd,r.exonCount,r.exonStarts,r.exonEnds,r.score,r.name2,r.cdsStartStat,r.cdsEndStat,r.exonFrames FROM refGene as r, hgFixed.gbCdnaInfo as i, ucscToEnsembl as c WHERE r.name=i.acc AND c.ucsc = r.chrom ORDER BY r.bin;" > refGene
