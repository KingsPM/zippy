all: install bowtie genome

install:
	apt-get update
	apt-get install -y python-pip python2.7-dev sqlite ncurses-dev unzip git python-virtualenv
	pip install cython
	pip install primer3-py
	pip install pysam
	pip install intervaltree

bowtie:
	wget http://netix.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip && \
	unzip bowtie2-2.2.6-linux-x86_64.zip && \
	cd bowtie2-2.2.6 && \
	mv bowtie2* /usr/local/bin
	rm -rf bowtie2-2.2.6

genome: genome-download genome-index

genome-download:
	mkdir -p zippy_resources && cd zippy_resources && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.fasta && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.fasta.fai

genome-index:
	cd zippy_resources && \
	bowtie2-build human_g1k_v37.fasta human_g1k_v37.bowtie

flask:
	virtualenv flask
	flask/bin/pip install flask
	flask/bin/pip install flask-login
	flask/bin/pip install flask-openid
	flask/bin/pip install flask-mail
	flask/bin/pip install flask-sqlalchemy
	flask/bin/pip install sqlalchemy-migrate
	flask/bin/pip install flask-whooshalchemy
	flask/bin/pip install flask-wtf
	flask/bin/pip install flask-babel
	flask/bin/pip install guess_language
	flask/bin/pip install flipflop
	flask/bin/pip install coverage
	mkdir tmp
