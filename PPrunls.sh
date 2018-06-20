#!/bin/bash

#PoolParty v0.81
#PPrunls


BASEDIR=$(dirname "$0")

if [[ ${BASEDIR} != *poolparty* ]];then
	SLOC=$(readlink -f "$0")
	BASEDIR=${SLOC%/*}
fi
#pass local score analysis to R script

set -o pipefail

#Define arguments
while getopts ":o:i:x:s:u:h:" opt; do
  case $opt in
    o) output="$OPTARG"
    ;;
    i) input="$OPTARG"
    ;;
    x) tuningparam="$OPTARG"
    ;;
    s) scaff="$OPTARG"
    ;;
    u) forceuni="$OPTARG"
    ;; 
    h) help="$OPTARG"
    ;;  
    \?) echo "Invalid option -$OPTARG" >&6
    ;;
  esac
done

#help
if [ -z "$output" ] || [ -z "$input" ] ; then
	printf "\nPPrunls : Uses local score R script to calculate significant regions from p-values (FET,CMH,FLK) \n \n"
	printf " Requires PPanalyze output format: 4 columns (Chr,pos,SnpID,p) \n "
	printf "Please read Accounting for linkage disequilibrium in genome scans forselection without individual genotypes: The local score approach (Fariello et al. 2017) \n\n "
	printf "Usage:   PPrunls -i [input] -o [output name] -s [scaffold prefix] -x [tuning parameter] \n"
	printf "\n Argument                             Description                                                     Default\n \n"
	printf " - i   [input]                        : Input file with four column format (CHR,POS,SNP,P)            [REQUIRED] \n"
	printf " - o   [output name]                  : Prefix for output name to define the analysis                 [REQUIRED] \n" 
	printf " - s   [scaffold prefix]              : If scaffolds present, prefix that identifies them             [scaff] \n"
	printf " - x   [tuning parameter]             : Override tuning parameter for local score (number)            [NULL] \n" 
	printf " - u   [force uniformity]             : Force p-value distribution to be recognized as uniform        [NULL] \n \n"
	exit
fi

[ -z "$scaff" ] && scaff=scaff
[ -z "$tuningparam" ] && tuningparam=NULL
[ -z "$forceuni" ] && forceuni=NULL


Rscript $BASEDIR/rscripts/r_localscore.R $input $output $scaff $tuningparam $forceuni
