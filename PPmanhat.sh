#!/bin/bash

#PoolParty v0.81
#PPmanhat


BASEDIR=$(dirname "$0")

if [[ ${BASEDIR} != *poolparty* ]];then
	SLOC=$(readlink -f "$0")
	BASEDIR=${SLOC%/*}
fi

set -o pipefail
set -e


#pass plotting options into R

#file to plot i
#chromosome ranges c
#log10 transformation l
#name o
#analysis type a


#Define arguments
while getopts ":o:i:c:l:a:s:r:C:R:1:2:L:p:t:m:z:h:" opt; do
  case $opt in
    o) output="$OPTARG"
    ;;
    i) input="$OPTARG"
    ;;
    c) ranges="$OPTARG"
    ;;
    l) logtrans="$OPTARG"
    ;;
    a) analtype="$OPTARG"
    ;;
    s) scaff="$OPTARG"
    ;;
    r) scale="$OPTARG"
    ;;
    C) chromosome="$OPTARG"
    ;;
    R) chrrange="$OPTARG"
    ;;
    1) color1="$OPTARG"
    ;;
    2) color2="$OPTARG"
    ;;
    L) gline="$OPTARG"
    ;;
    p) makepdf="$OPTARG"
    ;;
    t) plottype="$OPTARG"
    ;;
    m) minval="$OPTARG"
    ;;
    z) ztrans="$OPTARG"
    ;;
    h) help="$OPTARG"
    ;;  
    \?) echo "Invalid option -$OPTARG" >&17
    ;;
  esac
done

#help
if [ -z "$output" ] || [ -z "$input" ] ; then
	printf "\nPPManhat : Uses qqman to make Manhattan plots (.pdf, .png) of results from PPanalyze (FST,FET,SFST) \n \n"
	printf "Usage:   PPmanhat -i [input] -o [output name] -c [chromosome range file] -l [-log 10 trans] \n"
	printf "                  -a [analysis type] -s [scaffold heading] -r [scale by value] \n"
        printf "                  -C [plot chromosome] -R [plot range within chromosome] -1 [color1] -2 [color2] \n"
	printf " 	          -L [threshold line] -t [plot type] -p [make pdf too] -z [z transform values]"
	printf "\n Argument                             Description                                                     Default\n \n"
	printf " - i   [input]                        : Input file with four column format (CHR,POS,SNP,FST/P)        [REQUIRED] \n"
	printf " - o   [output name]                  : Prefix used to define the analysis                            [REQUIRED] \n" 
	printf " - c   [chromosome range file]        : File with chr names and starting/ending positions (CHRbp.txt) [NULL] \n"
	printf " - l   [-log 10 trans]                : Perform -log10 transform on column4                           [FALSE] \n" 
	printf " - a   [analysis type]                : Descriptor for Y axis (i.e., FST,-log10)                      [Differentiation] \n"
	printf " - s   [scaffold heading]             : String identifier for scaffolds, if present                   [scaffold] \n"
	printf " - r   [scale by value]               : If TRUE, Y range is scaled. If FALSE, Y range is 0-1          [TRUE] \n"
	printf " - C   [plot chromosome]              : Plot a single chromosome (integer)                            [NULL] \n"
	printf " - R   [plot range within chromosome] : Plot bp range within -C (two integers, comma-separated)       [NULL] \n"
	printf " - 1   [color1]                       : Color of odd chromosomes (R base colors)                      [darkred] \n"
	printf " - 2   [color2]                       : Color of even chromosomes (R base colors)                     [darkblue] \n"
	printf " - L   [threshold]                    : Draw threshold line at this value (integer)                   [NULL] \n"
	printf " - t   [plot type]                    : How to plot points (point,line,bar)                           [point] \n"
	printf " - m   [min value]                    : Don't plot SNPs less than this value (float/integer)          [0] \n"
	printf " - z   [z trans]                      : Z transform value                                             [FALSE] \n"
	printf " - p   [make pdf]                     : Make a pdf output in addition to png output                   [NULL] \n \n"
	exit
fi

[ -z "$ranges" ] && ranges=NULL
[ -z "$logtrans" ] && logtrans=NULL
[ -z "$analtype" ] && analtype=Differentiation
[ -z "$outdir" ] && outdir=$PWD
[ -z "$scaff" ] && scaff="scaffold"
[ -z "$scale" ] && scale=NULL
[ -z "$chromosome" ] && chromosome=NULL
[ -z "$chrrange" ] && chrrange=NULL
[ -z "$color1" ] && color1=darkred
[ -z "$color2" ] && color2=darkblue
[ -z "$gline" ] && gline=F
[ -z "$makepdf" ] && makepdf=NULL
[ -z "$plottype" ] && plottype=point
[ -z "$minval" ] && minval=0
[ -z "$ztrans" ] && ztrans=F


	Rscript $BASEDIR/rscripts/r_plotter.R $input $output $outdir $analtype $logtrans $ranges $scaff $scale $chromosome $chrrange $color1 $color2 $gline $makepdf $plottype $minval $ztrans

