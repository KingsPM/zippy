# installation makefile

all: install bowtie genome

install:
	locale-gen en_GB.UTF-8
	apt-get update
	apt-get install -y python-pip python2.7-dev sqlite3 ncurses-dev unzip git python-virtualenv htop
	apt-get install -y libxslt-dev libxml2-dev libffi-dev
	apt-get install -y redis-server
	apt-get install -y build-essential libjpeg-dev libfreetype6-dev python-dev python-imaging
	# create python venv
	virtualenv zippyvenv
	zippyvenv/bin/pip install reportlab
	zippyvenv/bin/pip install cython
	zippyvenv/bin/pip install primer3-py
	zippyvenv/bin/pip install pysam
	zippyvenv/bin/pip install intervaltree
	zippyvenv/bin/pip install flask
	zippyvenv/bin/pip install flask-login
	zippyvenv/bin/pip install flask-openid
	zippyvenv/bin/pip install flask-mail
	zippyvenv/bin/pip install flask-sqlalchemy
	zippyvenv/bin/pip install sqlalchemy-migrate
	zippyvenv/bin/pip install flask-whooshalchemy
	zippyvenv/bin/pip install flask-wtf
	zippyvenv/bin/pip install flask-babel
	zippyvenv/bin/pip install guess_language
	zippyvenv/bin/pip install flipflop
	zippyvenv/bin/pip install coverage
	zippyvenv/bin/pip install celery
	zippyvenv/bin/pip install redis
	mkdir tmp
	echo "activate python virtual environment with 'source zippyvenv/bin/activate'"

bowtie:
	wget -c http://netix.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip && \
	unzip bowtie2-2.2.6-linux-x86_64.zip && \
	cd bowtie2-2.2.6 && \
	mv bowtie2* /usr/local/bin
	rm -rf bowtie2-2.2.6

genome: genome-download genome-index

genome-download:
	mkdir -p zippy_resources && cd zippy_resources && \
	wget ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz && \
	wget ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai && \
	gunzip human_g1k_v37.fasta.gz

genome-index:
	cd zippy_resources && \
	bowtie2-build human_g1k_v37.fasta human_g1k_v37.bowtie
