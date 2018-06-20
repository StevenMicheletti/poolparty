#!/bin/bash

#PoolParty v0.81
#PPrunflk


BASEDIR=$(dirname "$0")

if [[ ${BASEDIR} != *poolparty* ]];then
	SLOC=$(readlink -f "$0")
	BASEDIR=${SLOC%/*}
fi

set -o pipefail

#pass an analyzed allele frequency table to FLK R script

#Define arguments
while getopts ":o:i:g:h:" opt; do
  case $opt in
    o) output="$OPTARG"
    ;;
    i) input="$OPTARG"
    ;;
    g) outgroup="$OPTARG"
    ;;
    h) help="$OPTARG"
    ;;  
    \?) echo "Invalid option -$OPTARG" >&4
    ;;
  esac
done


#help
if [ -z "$output" ] || [ -z "$input" ] ; then
	printf "\nPPrunflk : Runs FLK script using an allele frequency table (.fz) \n \n"
	printf "Usage:   PPrunflk -i [input] -o [output name] -g [outgroup] ] \n"
	printf "\n Argument                             Description                                                     Default\n \n"
	printf " - i   [input]                        : Input file in PP .fz format                                   [REQUIRED] \n"
	printf " - o   [output name]                  : Prefix for output name                                        [REQUIRED] \n" 
	printf " - g   [outgroup]                     : Specify an outgroup if desired (population integer)           [FALSE] \n"

	exit
fi

[ -z "$outgroup" ] && outgroup=NULL


Rscript $BASEDIR/rscripts/r_FLKu.R $input $output $outgroup


if [[ ! -s ${output}_flk_results.txt ]] ; then
	echo "ERROR: FLK failed, ensure you have > 2 populations"
	exit
fi

awk '{print $1,$2,NR, $7}'  ${output}_flk_results.txt | awk '{if (NR!=1) {print}}'  > ${output}.flk


echo "ALERT: Done"
