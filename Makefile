

install:
	apt-get update
	apt-get install -y python-pip python2.7-dev sqlite ncurses-dev unzip
	pip install cython
	pip install primer3-py
	pip install pysam
	pip install intervaltree


bowtie:
	wget http://netix.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip && \
	unzip bowtie2-2.2.6-linux-x86_64.zip && \
	cd bowtie2-2.2.6-linux-x86_64 && \
	make && make install 
# install BLAT

# get reference Genome
