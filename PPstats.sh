#!/bin/bash

#PoolParty v0.81
#PPstats

source $1
BASEDIR=$(dirname "$0")

if [[ ${BASEDIR} != *poolparty* ]];then
	SLOC=$(readlink -f "$0")
	BASEDIR=${SLOC%/*}
fi

set -o pipefail

if  ( [[ $(echo $1)  = "" ]] )  ; then
	echo "ERROR: You must provide a config file after the PPstats command"
	echo "for example, ./PPstats poolparty_stats.config"
	exit
fi
	echo "ALERT: $1 has been specified as the configuration file"

#Create output file
#This heading is intended for R to force 14 columns 

#Make the output directory if it does not exist 
	if [ ! -d "$OUTDIR/" ]; then
		mkdir $OUTDIR/
	fi

	if [[ -f ${OUTDIR}/${OUTFILE} ]] ; then
		echo "ERROR: ${OUTDIR}/${OUTFILE} already exists, cannot overwrite "
		exit
	fi

#Check for presence of MPILEUP

	if [[ ! -f ${MPILEUP} ]] ; then
   	 	echo "ERROR: mpileup file doesn't exist! Aborting."
    		exit
	fi

##SUBSET CHECK/OR RUN COMPLETE

#Ensure parameters are correct
#Ensure poolparty directory is intact
		if [ ! -d "$BASEDIR/rscripts" ] ; then
			echo "ERROR: Poolparty directory has been tampered with"
			exit
		fi

if [[ ! -f ${FAI} ]] ; then
    echo "ERROR: Genome .fai file doesn't exist! Aborting."
    exit
fi

if [[ ${FAI} =~ ".fai" ]];
    then
    :
else
    echo "ERROR: fai file does not have the correct extension. Genome fai file must end in .fai"
    exit
fi


if [[ ${MPILEUP} =~ "_stats" ]]; then
    :
else
    echo "ERROR: Mpileup should be the statistics version of the mpileup file which contains '_stats' in the filename "
    exit
fi



if [ ! -z "$SCAFP" ]; then
	if ! grep -q $SCAFP $FAI;
 	then
   		echo "ERROR: Scaffold prefix is wrong. Check fai file for scaffold prefix. Aborting"
 		exit
	fi
fi


if [ -z "$SCAFP" ]; then
     echo "ALERT: No scaffolds are identified and thus all scaffolds will be used according to the .fai file"
     SCAFP="ZcXYxxtxF3"
fi


if [[ $MINCOV =~ ^[\-0-9]+$ ]] && (( $MINCOV < 1)); then
    echo "ERROR: Minimum coverage must be a positive integer. Aborting."
 	exit
fi

if [[ $MAXCOV =~ ^[\-0-9]+$ ]] && (( $MINCOV < 1)); then
    echo "ERROR: Minimum coverage must be a positive integer. Aborting"
 	exit
fi

#Check for piping ability
if  [[ $(command -v mkfifo)  = "" ]] ; then
	echo "WARNING: piping 'mkfifo' not detected on system or available on drive. This may cause issues in downstream analyses"
	echo "Edit PPstats.sh and redirect mkfifo to another drive"
fi



#Check for R
		if  [[ $(command -v Rscript)  = "" ]] ; then
			echo "ERROR: R not detected on system. R is required for plotting "
			echo "R should initiate when 'Rscript' is typed in terminal"
			exit
		fi

#Load R to check for dependencies 
	if  [[ -f $OUTDIR/R_ERROR.txt ]] ; then
		rm $OUTDIR/R_ERROR.txt 
	fi		

	Rscript $BASEDIR/rscripts/r_stats_check.R $OUTDIR

	if  [[ -f $OUTDIR/R_ERROR.txt ]] ; then
		echo "ERROR: R dependency check failed, install dependencies manually"
		exit
	fi

echo "CONF: Parameter check passed. moving on..."

echo "V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14" > ${OUTDIR}/${OUTFILE}

#Check complete

#-Genome size determination (uses indexed file [.fai])-#
	#Anchored length
	AL=$(grep -v "^${SCAFP}" ${FAI} | awk '{ sum+=$2} END {printf "%.f\n", sum}') 
	#Full length
	FL=$(awk '{ sum+=$2} END {printf "%.f\n", sum}' ${FAI})
	#Scaffold length
	SL="$(($FL - $AL))"
	#Chromosomes
	CHRNUM=$(grep -v "^${SCAFP}" ${FAI} | cut -f 1  | uniq)
	CHRNUM2=$(grep -v "^${SCAFP}" ${FAI} | wc -l)

#Summary information on Genome
	echo "SUMMARY Full genome length is $FL bp" >> ${OUTDIR}/${OUTFILE}
	echo "SUMMARY Anchored genome length is $AL bp" >> ${OUTDIR}/${OUTFILE}
	echo "SUMMARY Scaffold length is $SL bp" >> ${OUTDIR}/${OUTFILE}
	echo "SCAFFOLD $SL" >> ${OUTDIR}/${OUTFILE}
	PROP=$(echo "$AL $FL" | awk '{printf "%.5f \n", $1/$2}')
	echo "SUMMARY Proportion of assembly that is anchored is $PROP" >> ${OUTDIR}/${OUTFILE}
	echo "SUMMARY There are $CHRNUM2 chromosomes" >> ${OUTDIR}/${OUTFILE}
	
	if (( ${CHRNUM2} > 100 )); then
	echo "WARNING more than 100 chromosomes, skipping chromosome plotting"
	fi

	echo "CHRNUM Chromosome bps" >> ${OUTDIR}/${OUTFILE}
	declare -a carr=($CHRNUM)
	for ((j=0;j<${#carr[@]};++j)); do
		CLENGTH=$(awk '$1 == "'${carr[j]}'" ' ${FAI} |  cut -f 1,2)
		echo "CHROMOSOME $CLENGTH" >> ${OUTDIR}/${OUTFILE}
	done

#-Column determination-#
	#Get Columns of mpileup and determine population number
		declare -i NCOL=$(awk '{ print NF; exit }' ${MPILEUP})
		declare -i NPOPS="$(($NCOL - 3))"
		POPSEQ=$(seq 4 1 $NCOL)
		#Select columns to run analyses on 
		declare -a arr2=($POPSEQ)
		POPNUM=$(seq 1 $NPOPS)

#-Print population summary-#
	echo "SUMMARY There are $NPOPS populations in the mpileup" >> ${OUTDIR}/${OUTFILE}
	echo "INFO Populations $(echo $POPNUM) correspond to columns $(echo $POPSEQ) in the mpileup" >> ${OUTDIR}/${OUTFILE}

#-Print population summary- done-#
	#Use sed to build an awk command that filters all populations at once
	ONE=$(echo $POPSEQ | sed -r 's/([^ ]+)/$\1/g')
	TWO=$(echo $ONE | awk '$NF=$NF " >= '$MINCOV' && "' OFS=" >= $MINCOV && ")
	THREE=$(echo $ONE | awk '$NF=$NF " <= '$MAXCOV'"' OFS=" <= $MAXCOV && ")
	FOUR=$(echo $TWO$THREE)
	
foo() {
		echo "STAT_TYPE Proportion of genome covered by all populations between $MINCOV X and $MAXCOV X" >> ${OUTDIR}/${OUTFILE}
			NLINE=$(awk  "$FOUR"  ${MPILEUP}  | wc -l)
			echo "COMB_TOT_BP $NLINE bp covered sufficiently by all libraries" >> ${OUTDIR}/${OUTFILE}
			PROP=$(echo "$NLINE $FL" | awk '{printf "%.5f \n", $1/$2}')
			echo "COMB_TOT_PROP $PROP of genome covered sufficiently by all libraries" >> ${OUTDIR}/${OUTFILE}
			echo "CONF 1/8 Complete, summary acquired and proportion of genomic overlap detected" 
}

	foo &

#HEADING  DONE
#SPAMMING STATS TO OUTPUT

##############
#MEAN COVERED
##############
	echo "STAT_TYPE Mean coverage of each population" >> ${OUTDIR}/${OUTFILE}
	#First column = pop 1

	
	open_sem(){
    		mkfifo pipe-$$
    		exec 3<>pipe-$$
  		  rm pipe-$$
   		 local i=$1
   		 for((;i>0;i--)); do
      			  printf %s 000 >&3
    		done
	}
	run_with_lock(){
   	 local x
   	 read -u 3 -n 3 x && ((0==x)) || exit $x
   	 (
    	"$@" 
   	 printf '%.3d' $? >&3
    	)&
	}

#Filtering bams; removing junk. Parallel.  
task() {

			ITER="$((${arr2[$i]} - 3))"
			DOLLA=$
			NEW="$DOLLA${arr2[i]}"
			MCOV=$(awk '{ sum += '$NEW'; n++ } END { if (n > 0) print sum / n; } ' ${MPILEUP})  
			MDEV=$(awk '{sum+='$NEW'; sumsq+='$NEW'*'$NEW'} END {print sqrt(sumsq/NR - (sum/NR)^2)}' ${MPILEUP})
			echo "TMA ${ITER} $MCOV" >> ${OUTDIR}/${OUTFILE}
			echo "TMS ${ITER} $MDEV" >> ${OUTDIR}/${OUTFILE}
			
 }

		N=$THREADZ
		open_sem $N
			for ((i=0;i<${#arr2[@]};++i)); do
				run_with_lock task $i
		done ; wait

echo "CONF 2/8 Complete, total mean coverage calculated"

############################
#MEAN COVERED AFTER FILTERS
#############################

	echo "STAT_TYPE Mean coverage of each populations between $MINCOV X and $MAXCOV X" >> ${OUTDIR}/${OUTFILE}


 task2() {
			ITER="$((${arr2[$i]} - 3))"
			#echo "Population $ITER"
			DOLLA=$
			NEW="$DOLLA${arr2[i]}"
			MCOV=$(awk ' '$NEW' >= '$MINCOV' && '$NEW' <= '$MAXCOV' ' ${MPILEUP} | awk '{ sum += '$NEW'; n++ } END { if (n > 0) print sum / n; } ')  
			MDEV=$(awk ' '$NEW' >= '$MINCOV' && '$NEW' <= '$MAXCOV' ' ${MPILEUP} | awk '{sum+='$NEW'; sumsq+='$NEW'*'$NEW'} END {print sqrt(sumsq/NR - (sum/NR)^2)}')
			echo "FMA ${ITER} $MCOV" >> ${OUTDIR}/${OUTFILE}
			echo "FMS ${ITER} $MDEV" >> ${OUTDIR}/${OUTFILE}
 }

		N=$THREADZ
		open_sem $N
			for ((i=0;i<${#arr2[@]};++i)); do
				run_with_lock task2 $i
		done ; wait

echo "CONF 3/8 Complete, filtered mean coverage calculated"

##############
#PROP COVERED
##############

	echo "STAT_TYPE Proportion of genome covered between $MINCOV X and $MAXCOV X " >> ${OUTDIR}/${OUTFILE}

 task3() {	
			ITER="$((${arr2[$i]} - 3))"
			DOLLA=$
			NEW="$DOLLA${arr2[i]}"
			NLINE=$(awk ' '$NEW' >= '$MINCOV' && '$NEW' <= '$MAXCOV'  ' ${MPILEUP}  | wc -l)
			echo "FBC ${ITER} $NLINE" >> ${OUTDIR}/${OUTFILE}
			PROP=$(echo "$NLINE $FL" | awk '{printf "%.5f \n", $1/$2}')
			echo "FPC ${ITER} $PROP" >> ${OUTDIR}/${OUTFILE}
 }

		N=$THREADZ
		open_sem $N
			for ((i=0;i<${#arr2[@]};++i)); do
				run_with_lock task3 $i
		done ; wait

echo "CONF 4/8 Complete, proportion of ref genome calculated"

##############################
#Proportion Scaffold  Covered
##############################

CHRNUM=$(grep "^${SCAFP}" ${FAI}| cut -f 1  | uniq)
	echo "STAT_TYPE Scaffold proportion of alignments between $MINCOV and $MAXCOV" >> ${OUTDIR}/${OUTFILE}
	
 task4()	{
			ITER="$((${arr2[$j]} - 3))"
			#echo "Population $ITER for $MINCOV to $MAXCOV"
			DOLLA=$
			COL=${arr2[j]}
			NEW="$DOLLA$COL"
			NLINE=$(grep "^${SCAFP}" ${MPILEUP} | awk '  '$NEW' >= '$MINCOV' && '$NEW' <= '$MAXCOV' ' | wc -l)
			echo "SCAFL ${ITER} scaff $NLINE" >> ${OUTDIR}/${OUTFILE}
			PROP=$(echo "$NLINE $SL" | awk '{printf "%.5f \n", $1/$2}')
			echo "SCAFP ${ITER} scaff $PROP" >> ${OUTDIR}/${OUTFILE}
 }

		N=$THREADZ
		open_sem $N
			for ((j=0;j<${#arr2[@]};++j)); do
				run_with_lock task4 $j
		done ; wait

echo "CONF 5/8 Complete, scaffold proportion calculated"

##############################
#Proportion Anchored Covered
##############################

	echo "STAT_TYPE Anchored proportion of alignments between $MINCOV and $MAXCOV" >> ${OUTDIR}/${OUTFILE}
	
 task5()	{
			ITER="$((${arr2[$j]} - 3))"
			#echo "Population $ITER for $MINCOV to $MAXCOV"
			DOLLA=$
			COL=${arr2[j]}
			NEW="$DOLLA$COL"
			NLINE=$(grep -v "^${SCAFP}" ${MPILEUP} | awk ' '$NEW' >= '$MINCOV' && '$NEW' <= '$MAXCOV' ' | wc -l)
			echo "ANCL ${ITER} anch $NLINE" >> ${OUTDIR}/${OUTFILE}
			PROP=$(echo "$NLINE $AL" | awk '{printf "%.5f \n", $1/$2}')
			echo "ANCP ${ITER} anch $PROP" >> ${OUTDIR}/${OUTFILE}
 }

		N=$THREADZ
		open_sem $N
			for ((j=0;j<${#arr2[@]};++j)); do
				run_with_lock task5 $j
		done ; wait

echo "CONF 6/8 Complete, anchored proportion calculated"

###############################################
#Coverage at threshold intervals- outer loop #
###############################################

       echo "STAT_TYPE Proportion of full genome covered at different depths of coverage" >> ${OUTDIR}/${OUTFILE}
	for j in {1,5,10,15,20,25,30,40,50,75,100,150} ; do 
	

	#IND POP COVER THRESHOLD - FULL 
	#First column = pop 1
	task6() {
		
			ITER="$((${arr2[$i]} - 3))"
			DOLLA=$
			NEW="$DOLLA${arr2[i]}"
			NLINE=$(awk ' '$NEW' >= '$j' ' ${MPILEUP} | wc -l)
			echo "DOC ${j} ${ITER} $NLINE" >> ${OUTDIR}/${OUTFILE}
			PROP=$(echo "$NLINE $FL" | awk '{printf "%.5f \n", $1/$2}')
			echo "DOP ${j} ${ITER} $PROP" >> ${OUTDIR}/${OUTFILE}

		 }

 		N=$THREADZ
		open_sem $N
			for ((i=0;i<${#arr2[@]};++i)); do
				run_with_lock task6 $i
		done 

done ; wait

echo "CONF 7/8 Complete, interval proportion calculated"	

#####################################
#Chromosome-by-chromosome analysis#
#####################################

if (( ${CHRNUM2} < 10000 )); then

	CHRNUM=$(grep -v "^${SCAFP}" ${FAI}| cut -f 1  |uniq)
		echo "STAT_TYPE Chromosome-by-chromosome proportion of alignments between $MINCOV and $MAXCOV" >> ${OUTDIR}/${OUTFILE}
		declare -a carr=($CHRNUM)
		for ((j=0;j<${#arr2[@]};++j)); do
				ITER="$((${arr2[$j]} - 3))"
				DOLLA=$
				COL=${arr2[j]}
				NEW="$DOLLA$COL"
task7() {		

			CLENGTH=$(awk '$1 == "'${carr[i]}'" ' ${FAI} |  cut -f 2)
			NLINE=$(awk '$1 == "'${carr[i]}'" && '$NEW' >= '$MINCOV' && '$NEW' <= '$MAXCOV' ' ${MPILEUP}  | wc -l)
			echo "CHRC ${ITER} ${carr[i]} $NLINE" >> ${OUTDIR}/${OUTFILE}
			PROP=$(echo "$NLINE $CLENGTH" | awk '{printf "%.5f \n", $1/$2}')
			echo "CHRP ${ITER} ${carr[i]} $PROP" >> ${OUTDIR}/${OUTFILE}
}
 		N=$THREADZ
		open_sem $N
			for ((i=0;i<${#carr[@]};++i)); do
				run_with_lock task7 $i
		done 
	done ; wait

echo "CONF 8/8 Complete, chromosome proportion calculated"

fi 
echo "CONF BASH portion of run is complete"
echo "CONF Loading R to plot results and create tables"
	rin=${OUTDIR}/${OUTFILE}
	rout=${OUTDIR}
		Rscript ${BASEDIR}/rscripts/r_plotstats.R $rin $rout

echo "CONF all files have been written to ${OUTDIR}"