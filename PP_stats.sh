#!/bin/bash
source $PWD/PP_stats.conf

echo V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14
echo V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14

#Check for presence of MPILEUP

	if [[ ! -f ${MPILEUP} ]] ; then
   	 	echo "mpileup file doesn't exist! Aborting."
    		exit
	fi

	if [[ ${MPILEUP} != *"_stats.mp"* ]];then
	echo "ALERT mpileup is not reduced. Filtering now"
   	 #MPILEUP has not been filtered, beginning subset process
		#Get MPILEUP new output name
			MPO=$(echo ${MPILEUP} | awk -F'[.]' '{print $1}')
			MPN=${MPO}_stats.mp
		#Get Columns of mpileup and determine population number
			declare -i NCOL=$(awk '{ print NF; exit }' ${MPILEUP})
			#population number (3 columns per population)
			declare -i NPOPS=($NCOL-3)/3
		# Start at 3, increment by NPOPS until NCOL, i.e, get column of each pop
			POPSEQ=$(seq 4 3 $NCOL)
			POPNUM=$(seq 1 $NPOPS)
		#Select columns to run analyses on 
			declare -a arr2=($POPSEQ)
		#Subset
			ONE=$(echo $POPSEQ | sed -r 's/([^ ]+)/$\1/g')
			TWO=$(echo $ONE | awk '$NF=$NF "}" ' OFS=",")
			THREE='{print $1,$2,'
			FOUR=$(echo ${THREE}${TWO})	
			awk "$FOUR" ${MPILEUP} > ${MPN}
			declare -i NCOL2=$(awk '{ print NF; exit }' ${MPN})
			POPSEQ=$(seq 3 1 $NCOL2)
		#Select columns to run analyses on 
			declare -a arr2=($POPSEQ)
			NPOPS=`expr $NCOL2 - 2`
			MPILEUP=${MPN}
	fi


##SUBSET CHECK/OR RUN COMPLETE

#Ensure parameters are correct
if [[ ! -f ${FAI} ]] ; then
    echo ".fai file doesn't exist! Aborting."
    exit
fi

if ! grep -q "$SCAFP" $FAI;
 then
    echo "Scaffold prefix is wrong. Check fai file for scaffold prefix. Aborting"
 	exit
fi

if [[ $MINCOV =~ ^[\-0-9]+$ ]] && (( $MINCOV < 1)); then
    echo "Minimum coverage must be a positive integer. Aborting."
 	exit
fi

if [[ $MAXCOV =~ ^[\-0-9]+$ ]] && (( $MINCOV < 1)); then
    echo "Minimum coverage must be a positive integer. Aborting"
 	exit
fi

#Check complete

#-Genome size determination (uses indexed file [.fai])-#
	#Anchored length
	AL=$(grep -v "^${SCAFP}" ${FAI} | awk '{ sum+=$2} END {print sum}') 
	#Full length
	FL=$(awk '{ sum+=$2} END {print sum}' ${FAI})
	#Scaffold length
	SL=`expr $FL - $AL`
	#Chromosomes
	CHRNUM=$(grep -v "^${SCAFP}" ${FAI} | cut -f 1  | uniq)
	CHRNUM2=$(grep -v "^${SCAFP}" ${FAI} | wc -l)

#Summary information on Genome
	echo "SUMMARY Full genome length is $FL bp"
	echo "SUMMARY Anchored genome length is $AL bp"
	echo "SUMMARY Scaffold length is $SL bp"
	echo "SCAFFOLD $SL"
	PROP=$(echo "$AL $FL" | awk '{printf "%.5f \n", $1/$2}')
	echo "SUMMARY Proportion of assembly that is anchored is $PROP"
	echo "SUMMARY There are $CHRNUM2 chromosomes"
	
	if (( ${CHRNUM2} > 100 )); then
	echo "WARNING more than 100 chromosomes, skipping chromosome analyses"
	fi

	echo "CHRNUM Chromosome bps"
	declare -a carr=($CHRNUM)
	for ((j=0;j<${#carr[@]};++j)); do
		CLENGTH=$(awk '$1 == "'${carr[j]}'" ' ${FAI} |  cut -f 1,2)
		echo "CHROMOSOME $CLENGTH"
	done

#-Column determinination-#
	#Get Columns of mpileup and determine population number
		declare -i NCOL=$(awk '{ print NF; exit }' ${MPILEUP})
		declare -i NPOPS=`expr $NCOL - 2`
		POPSEQ=$(seq 3 1 $NCOL)
		#Select columns to run analyses on 
		declare -a arr2=($POPSEQ)

#-Print population summary-#
	echo "SUMMARY There are $NPOPS populations in the mpileup"
	echo INFO Populations $POPNUM correspond to columns $POPSEQ in the mpileup

#-Print population summary- done-#
	#Use sed to build an awk command that filters all populations at once
	ONE=$(echo $POPSEQ | sed -r 's/([^ ]+)/$\1/g')
	TWO=$(echo $ONE | awk '$NF=$NF " >= '$MINCOV' && "' OFS=" >= $MINCOV && ")
	THREE=$(echo $ONE | awk '$NF=$NF " <= '$MAXCOV'"' OFS=" <= $MAXCOV && ")
	FOUR=$(echo $TWO$THREE)
	
		echo "STAT_TYPE Proportion of genome covered by all populations between $MINCOV X and $MAXCOV X"
			NLINE=$(awk  "$FOUR"  ${MPILEUP}  | wc -l)
			echo "COMB_TOT_BP $NLINE bp covered sufficiently by all libraries"
			PROP=$(echo "$NLINE $FL" | awk '{printf "%.5f \n", $1/$2}')
			echo "COMB_TOT_PROP $PROP of genome covered sufficiently by all libraries"


echo "CONF 1/8 Complete"

#HEADING  DONE
#SPAMMING STATS TO OUTPUT

##############
#MEAN COVERED
##############
	echo "STAT_TYPE Mean coverage of each population"
	#First column = pop 1


 foo() {

			ITER=`expr ${arr2[$i]} - 2`
			DOLLA=$
			NEW="$DOLLA${arr2[i]}"
			MCOV=$(awk '{ sum += '$NEW'; n++ } END { if (n > 0) print sum / n; } ' ${MPILEUP})  
			MDEV=$(awk '{sum+='$NEW'; sumsq+='$NEW'*'$NEW'} END {print sqrt(sumsq/NR - (sum/NR)**2)}' ${MPILEUP})
			echo "TMA ${ITER} $MCOV"
			echo "TMS ${ITER} $MDEV"
			
 }

 for ((i=0;i<${#arr2[@]};++i)); do foo $i & done; wait

echo "CONF 2/8 Complete"

############################
#MEAN COVERED AFTER FILTERS
#############################

	echo "STAT_TYPE Mean coverage of each populations between $MINCOV X and $MAXCOV X"


 fooc() {
			ITER=`expr ${arr2[$i]} - 2`
			#echo "Population $ITER"
			DOLLA=$
			NEW="$DOLLA${arr2[i]}"
			MCOV=$(awk ' '$NEW' >= '$MINCOV' && '$NEW' <= '$MAXCOV' ' ${MPILEUP} | awk '{ sum += '$NEW'; n++ } END { if (n > 0) print sum / n; } ')  
			MDEV=$(awk ' '$NEW' >= '$MINCOV' && '$NEW' <= '$MAXCOV' ' ${MPILEUP} | awk '{sum+='$NEW'; sumsq+='$NEW'*'$NEW'} END {print sqrt(sumsq/NR - (sum/NR)**2)}')
			echo "FMA ${ITER} $MCOV"
			echo "FMS ${ITER} $MDEV"
 }

 for ((i=0;i<${#arr2[@]};++i)); do fooc $i & done; wait

echo "CONF 3/8 Complete"

##############
#PROP COVERED
##############

	echo "STAT_TYPE Proportion of genome covered between $MINCOV X and $MAXCOV X "

 fooa () {	ITER=`expr ${arr2[$i]} - 2`

			DOLLA=$
			NEW="$DOLLA${arr2[i]}"
			NLINE=$(awk ' '$NEW' >= '$MINCOV' && '$NEW' <= '$MAXCOV'  ' ${MPILEUP}  | wc -l)
			echo "FBC ${ITER} $NLINE"
			PROP=$(echo "$NLINE $FL" | awk '{printf "%.5f \n", $1/$2}')
			echo "FPC ${ITER} $PROP"
 }

 for ((i=0;i<${#arr2[@]};++i)); do fooa $i & done; wait

echo "CONF 4/8 Complete"

##############################
#Proportion Scaffold  Covered
##############################

CHRNUM=$(grep "^${SCAFP}" ${FAI}| cut -f 1  | uniq)
	echo "STAT_TYPE Scaffold proportion of alignments between $MINCOV and $MAXCOV"
	
 foog()	{
			ITER=`expr ${arr2[$j]} - 2`
			#echo "Population $ITER for $MINCOV to $MAXCOV"
			DOLLA=$
			COL=${arr2[j]}
			NEW="$DOLLA$COL"
			NLINE=$(grep "^${SCAFP}" ${MPILEUP} | awk '  '$NEW' >= '$MINCOV' && '$NEW' <= '$MAXCOV' ' | wc -l)
			echo "SCAFL ${ITER} scaff $NLINE"
			PROP=$(echo "$NLINE $SL" | awk '{printf "%.5f \n", $1/$2}')
			echo "SCAFP ${ITER} scaff $PROP"
 }

 for ((j=0;j<${#arr2[@]};++j)); do foog $j & done; wait

echo "CONF 5/8 Complete"



##############################
#Proportion Anchored Covered
##############################

	echo "STAT_TYPE Anchored proportion of alignments between $MINCOV and $MAXCOV"
	
 fooj()	{
			ITER=`expr ${arr2[$j]} - 2`
			#echo "Population $ITER for $MINCOV to $MAXCOV"
			DOLLA=$
			COL=${arr2[j]}
			NEW="$DOLLA$COL"
			NLINE=$(grep -v "^${SCAFP}" ${MPILEUP} | awk ' '$NEW' >= '$MINCOV' && '$NEW' <= '$MAXCOV' ' | wc -l)
			echo "ANCL ${ITER} anch $NLINE"
			PROP=$(echo "$NLINE $AL" | awk '{printf "%.5f \n", $1/$2}')
			echo "ANCP ${ITER} anch $PROP"
 }

 for ((j=0;j<${#arr2[@]};++j)); do fooj $j & done; wait

echo "CONF 6/8 Complete"	


###############################################
#Coverage at threshold intervals- outer loop #
###############################################

       echo "STAT_TYPE Proportion of full genome covered at different depths of coverage"
	for j in {1,5,10,15,20,25,30,40,50,75,100,150} ; do 
	

	#IND POP COVER THRESHOLD - FULL 
	#First column = pop 1
	food() {
		
			ITER=`expr ${arr2[$i]} - 2`
			DOLLA=$
			NEW="$DOLLA${arr2[i]}"
			NLINE=$(awk ' '$NEW' >= '$j' ' ${MPILEUP} | wc -l)
			echo "DOC ${j} ${ITER} $NLINE"
			PROP=$(echo "$NLINE $FL" | awk '{printf "%.5f \n", $1/$2}')
			echo "DOP ${j} ${ITER} $PROP"

		 }

 for ((i=0;i<${#arr2[@]};++i)); do food $i & done

done ; wait

echo "CONF 7/8 Complete"	

#####################################
#Chromosome-by-chromosome analysis#
#####################################

if (( ${CHRNUM2} < 100 )); then

CHRNUM=$(grep -v "^${SCAFP}" ${FAI}| cut -f 1  |uniq)
	echo "STAT_TYPE Chromosome-by-chromosome proportion of alignments between $MINCOV and $MAXCOV"
	declare -a carr=($CHRNUM)
	for ((j=0;j<${#arr2[@]};++j)); do
			ITER=`expr ${arr2[$j]} - 2`
			DOLLA=$
			COL=${arr2[j]}
			NEW="$DOLLA$COL"
foof() {		

			CLENGTH=$(awk '$1 == "'${carr[i]}'" ' ${FAI} |  cut -f 2)
			NLINE=$(awk '$1 == "'${carr[i]}'" && '$NEW' >= '$MINCOV' && '$NEW' <= '$MAXCOV' ' ${MPILEUP}  | wc -l)
			echo "CHRC ${ITER} ${carr[i]} $NLINE"
			PROP=$(echo "$NLINE $CLENGTH" | awk '{printf "%.5f \n", $1/$2}')
			echo "CHRP ${ITER} ${carr[i]} $PROP"
}
for ((i=0;i<${#carr[@]};++i)); do foof $i & done
	done ; wait

echo "CONF 8/8 Complete"	


fi


echo "CONF Full run is complete"	

