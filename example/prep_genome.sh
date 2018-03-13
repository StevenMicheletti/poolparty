#!/bin/bash

#prepares reference genome in FASTA format for BWA alignment

GENOME="PP_genome.fa"	

#!/bin/bash
# 1) Index fastA genome file with bwa
	
bwa index -a bwtsw ${GENOME}

# 2) Create fastA file index as well

samtools faidx ${GENOME}

# 3) Create sequence dictionary with picardtools (Java)

java -Xmx2g -jar /usr/local/bin/picard.jar CreateSequenceDictionary REFERENCE= ${GENOME} OUTPUT=${GENOME}.dict

