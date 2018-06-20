#bin/bash!

#PoolParty v0.81
#PPruncmh


BASEDIR=$(dirname "$0")

if [[ ${BASEDIR} != *poolparty* ]];then
	SLOC=$(readlink -f "$0")
	BASEDIR=${SLOC%/*}
fi



#Define arguments
while getopts ":P:i:o:p:h:m:x:" opt; do
  case $opt in
    P) pooldir="$OPTARG"
    ;;
    i) input="$OPTARG"
    ;;
    o) outname="$OPTARG"
    ;;
    p) pops="$OPTARG"
    ;;
    h) help="$OPTARG"
    ;;
    m) mincov="$OPTARG"
    ;;
    x) maxcov="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&7
    ;;
  esac
done


#help
if [ -z "$pooldir" ] || [ -z "$input" ] || [ -z "$outname" ] ; then
	printf "\nRun Popoolation2 CMH test on specified populations \n \n"
	printf "Usage: PPruncmh -i [input] -o [outname] -d [outdir] \n"
	printf "                -p [pops] -P [pooldir] \n \n"
	printf "Argument            Description                                                       Default  \n \n"
	printf " - i    [input]    : Input sync file                                                  [REQUIRED] \n"
	printf " - o    [outname]  : Output prefix name                                               [REQUIRED] \n" 
	printf " - p    [pops]     : pops to compare. i.e., 1-2,3-4  where 1&2 are from same locality [REQUIRED] \n" 
	printf " - P    [pooldir]  : Popoolation2 base directory                                      [REQUIRED] \n"
	printf " - m    [mincov]   : Minimum coverage override for unanalyzed sync file               [1] \n"
	printf " - x    [maxcov]   : Maximum coverage override for unanalyzed sync file               [10000] \n \n"
	echo "Check PPanalyze for sync population order"
	exit
fi

[ -z "$outdir" ] && outdir=$PWD
[ -z "$mincov" ] && mincov=1
[ -z "$maxcov" ] && maxcov=10000

echo "ALERT: Running CMH test at $(date)"

perl ${pooldir}/cmh-test.pl --input ${input} --output ${outdir}/${outname}_raw.cmh --min-count 1 --min-coverage ${mincov} --max-coverage ${maxcov} --population ${pops}

echo "ALERT: CMH test done at $(date)"
echo "ALERT: Formatting CMH results"

	paste -d' ' <(awk {'print $1,$2'} ${outdir}/${outname}_raw.cmh ) <(awk '{print NR}' ${outdir}/${outname}_raw.cmh) <(awk '{print $(NF-0)}' ${outdir}/${outname}_raw.cmh) > ${outdir}/${outname}.cmh
	rm ${outdir}/${outname}_raw.cmh.rin
	rm ${outdir}/${outname}_raw.cmh.rout
	rm ${outdir}/${outname}_raw.cmh.params
	
echo "ALERT: Done with CMH test"