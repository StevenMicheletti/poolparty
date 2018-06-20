#bin/bash!

#PoolParty v0.81
#PPsubset


BASEDIR=$(dirname "$0")

if [[ ${BASEDIR} != *poolparty* ]];then
	SLOC=$(readlink -f "$0")
	BASEDIR=${SLOC%/*}
fi

set -o pipefail
set -e

#Define arguments
while getopts ":d:m:x:o:P:n:t:a:p:h:" opt; do
  case $opt in
    d) output="$OPTARG"
    ;;
    m) mincov="$OPTARG"
    ;;
    x) maxcov="$OPTARG"
    ;;
    o) outname="$OPTARG"
    ;;
    P) propsnps="$OPTARG"
    ;;
    n) numsnps="$OPTARG"
    ;;
    t) analtype="$OPTARG"
    ;;
    a) analname="$OPTARG"
    ;;
    p) pops="$OPTARG"
    ;;
    h) help="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&11
    ;;
  esac
done


#help
if [ -z "$analtype" ] || [ -z "$analname" ] ; then
	printf "\nSubset PPanalyze results by various thresholds \n \n"
	printf "Usage: PPsubset -t [analtype] -a [analname] -o [output] -p [pops] \n"
	printf "                -m [mincov] -x [maxcov] -P [propsnps] -n [numsnps] \n \n"
	printf "Argument            Description                            Default  \n \n"
	printf " - t    [analtype]  : either FST, SFST, or FET             [REQUIRED] \n"
	printf " - a    [analname]  : prefix used to define the analysis   [REQUIRED] \n" 
	printf " - d    [directory] : location of analyzed files           [CURRENTDIR] \n"
	printf " - o    [outname]   : Id suffix of reclassed files         [rc] \n"
	printf " - m    [mincov]    : minimum coverage                     [1] \n" 
	printf " - x    [maxcov]    : maximum coverage                     [10000] \n"
	printf " - P    [propsnps]  : minimum proportion of SNPs in window [0] \n"
	printf " - n    [numsnps]   : minimum number of SNPs in window     [1] \n"
	printf " - p    [pops]      : reclassification of populations      [NULL] \n\n "

	echo "If > 2 populations, script requires analname_avoidcols.txt"
	echo " Pop classification uses same format as PPanalyze"
	
	exit
fi

[ -z "$mincov" ] && mincov=1
[ -z "$maxcov" ] && maxcov=10000
[ -z "$propsnps" ] && propsnps=0
[ -z "$numsnps" ] && numsnps=1
[ -z "$output" ] && output=$PWD
[ -z "$pops" ] && pops=NULL
[ -z "$outname" ] && outname=rc


	if [ ! -d "${output}/temp" ]; then
		mkdir ${output}/temp
	fi
	

#If populations are being reclassified, create new avoid columns file

if [[ "$pops" != "NULL" ]] ; then

	POPS2=$(echo $pops | sed 's/[,:]/ /g')
		declare -a arr2=($POPS2)
		declare -i NPOPS=${#arr2[@]}
	
	echo "ALERT: There are $NPOPS populations in the new pop classification"
	
	if [ -z "$(echo $pops | sed -n 's/\([,:]\)/\1/p')" ];
		then echo "ERROR: pops specification contains invalid characters"
		exit
	fi
	
	#Create matching sequence for column name adjustment
	pseq=$(printf '%0.s3 ' $(seq 1 $NPOPS))

	
	if  [[ -f ${output}/${analname}_avoidcols_${outname}.txt ]] ; then
		rm ${output}/${analname}_avoidcols_${outname}.txt
	fi	
		
	#Prepare expression for keeping columns
		echo $pops | tr , \\n > ${output}/${analname}_avoidcols_${outname}.txt
		
	echo "ALERT: Using pops $POPS2 from the subset analysis. These correspond to the pop order in the analyzed file"
	
else
	if  [[ -f ${output}/${analname}_avoidcols.txt ]] ; then
		cp  ${output}/${analname}_avoidcols.txt ${output}/${analname}_avoidcols_${outname}.txt 
	fi
fi
	
	
#reclassify FST based on coverage 

if [[ "$analtype" == "FST" ]]; then
	echo "ALERT: Subsetting FST file ${analname}_raw.fst"
	awk '$5 >= '$mincov' &&  $5  <= '$maxcov' && $4 >= '$propsnps' && $3 >= '$numsnps' ' $output/${analname}_raw.fst >  $output/temp/${analname}_reclassing${outname}.fst

				declare -i before=$(wc -l $output/${analname}_raw.fst |  cut -f1 -d' ')
				declare -i afterc=$(wc -l $output/temp/${analname}_reclassing${outname}.fst | cut -f1 -d' ')
				declare -i lossrc=$(($before-$afterc))
				echo "ALERT: $lossrc SNPs were removed from FST analysis due to new parameters"


	#Cut up FST file
	cut -d$'\t' -f 1-2 $output/temp/${analname}_reclassing${outname}.fst | gawk '$3=(FNR FS $3)' > $output/temp/${analname}_head${outname}.fst
	cut -d$'\t' -f 6- $output/temp/${analname}_reclassing${outname}.fst > $output/temp/${analname}_body${outname}.fst
	awk 'NR==1{print $0}' $output/temp/${analname}_body${outname}.fst | awk  '{gsub("=","\t",$0); print;}' |  awk '{ for (i=2;i<=NF;i+=2) $i="" } 1' | awk  '{gsub("  ","\t",$0); print;}' | awk  '{gsub(" ","",$0); print;}'   > $output/temp/${analname}_body_heading${outname}.fst
	awk '{gsub("=","\t",$0); print;}' $output/temp/${analname}_body${outname}.fst | awk '{ for (i=1;i<=NF;i+=2) $i="" } 1' | awk  '{gsub("  ","\t",$0); print;}' | awk  '{gsub(" ","",$0); print;}'  > $output/temp/${analname}_body2${outname}.fst
	cat $output/temp/${analname}_body_heading${outname}.fst $output/temp/${analname}_body2${outname}.fst > $output/temp/${analname}_prep${outname}.fst

	if  [[ -f ${output}/${analname}_avoidcols_${outname}.txt ]] ; then
		if [[ "$pops" != "NULL" ]] ; then
			COMBLIST2=$(awk '{print}' ORS=' ' ${output}/${analname}_avoidcols_${outname}.txt) 
			COMBLIST3=$(echo $COMBLIST2 | sed 's/[^ ][^ ]*/"&"/g')
				ZERO=$(echo $COMBLIST3 | sed -r 's/([^ ]+)/$i==\1/g') 
				TWO=$(echo $ZERO | sed -e 's/ /\  \|\|\ /g')
				ONE=$(echo " NR==1{for(i=1; i<=NF; i++) if (")
				THREE=$(echo ' ) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"} ')
				FOUR=$(echo ${ONE}${TWO}${THREE})
		else
			COMBLIST2=$(awk '{print}' ORS=' ' ${output}/${analname}_avoidcols_${outname}.txt) 
			COMBLIST3=$(echo $COMBLIST2 | sed 's/[^ ][^ ]*/"&"/g')
				ZERO=$(echo $COMBLIST3 | sed -r 's/([^ ]+)/$i!=\1/g') 
				TWO=$(echo $ZERO | sed -e 's/ /\  \&\&\ /g')
				ONE=$(echo " NR==1{for(i=1; i<=NF; i++) if (")
				THREE=$(echo ' ) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"} ')
				FOUR=$(echo ${ONE}${TWO}${THREE})
		fi
			#Use above expression to remove FST columns in similar populations
			awk "$FOUR" $output/temp/${analname}_prep${outname}.fst > $output/temp/${analname}_keep${outname}.fst
			NCOLZ=$(awk '{print NF}' $output/temp/${analname}_keep${outname}.fst | sort -nu | tail -n 1)
			echo $NCOLZ
			if [ "$NCOLZ" -gt 1 ]; then
				awk '{if (NR!=1) {print}}' $output/temp/${analname}_keep${outname}.fst | awk '{ s = 0; for (i = 1; i <= NF; i++) s += $i; print $1, (NF > 1) ? s / (NF - 0) : 0; }' | awk '{print $2}' > $output/temp/${analname}_keep2${outname}.fst
			else
				awk '{if (NR!=1) {print}}' $output/temp/${analname}_keep${outname}.fst > $output/temp/${analname}_keep2${outname}.fst
			fi
			#Create final fst file for this analysis.
			paste $output/temp/${analname}_head${outname}.fst $output/temp/${analname}_keep2${outname}.fst | awk '{gsub("\t","",$0); print;}'  > $output/temp/${analname}_fstNA${outname}.fst
	else 
		declare -i NCOLZ=$(awk '{print NF}' $output/temp/${analname}_prep${outname}.fst | sort -nu | tail -n 1)
		if [ "$NCOLZ" -gt 1 ]; then
			awk '{if (NR!=1) {print}}' $output/temp/${analname}_prep${outname}.fst | awk '{ s = 0; for (i = 1; i <= NF; i++) s += $i; print $1, (NF > 1) ? s / (NF - 0) : 0; }' | awk '{print $2}' > $output/temp/${analname}_keep2${outname}.fst
		else
			awk '{if (NR!=1) {print}}' $output/temp/${analname}_prep${outname}.fst > $output/temp/${analname}_keep2${outname}.fst
		fi
			#Create final fst file for this analysis.
			paste $output/temp/${analname}_head${outname}.fst $output/temp/${analname}_keep2${outname}.fst | awk '{gsub("\t","",$0); print;}'  > $output/temp/${analname}_fstNA${outname}.fst
	fi

				awk '!/na/' $output/temp/${analname}_fstNA${outname}.fst > $output/${analname}_${outname}.fst
				declare -i before=$(wc -l $output/${analname}_raw.fst |  cut -f1 -d' ')
				declare -i after=$(wc -l $output/${analname}_${outname}.fst | cut -f1 -d' ')
				declare -i loss=$(($before-$after))
				#echo "ALERT: $loss SNPs were removed from FST analysis due to Ns, indels, or uninformative comparisons"
				echo "ALERT: $after SNPs remain"
				rm $output/temp/*${analname}*.fst*
fi

if [[ "$analtype" == "SFST" ]]; then
	echo "ALERT: Subsetting FST file ${analname}_raw.Sfst"
	awk '$5 >= '$mincov' &&  $5  <= '$maxcov' && $4 >= '$propsnps' && $3 >= '$numsnps' ' $output/${analname}_raw.Sfst >  $output/temp/${analname}_reclassing${outname}.Sfst

				declare -i before=$(wc -l $output/${analname}_raw.Sfst |  cut -f1 -d' ')
				declare -i afterc=$(wc -l $output/temp/${analname}_reclassing${outname}.Sfst | cut -f1 -d' ')
				declare -i lossrc=$(($before-$afterc))
				echo "ALERT: $lossrc SNP positions were removed from SFST analysis due to new parameters"


	#Cut up SFST file
	cut -d$'\t' -f 1-2 $output/temp/${analname}_reclassing${outname}.Sfst | gawk '$3=(FNR FS $3)' > $output/temp/${analname}_head${outname}.Sfst
	cut -d$'\t' -f 6- $output/temp/${analname}_reclassing${outname}.Sfst > $output/temp/${analname}_body${outname}.Sfst
	awk 'NR==1{print $0}' $output/temp/${analname}_body${outname}.Sfst | awk  '{gsub("=","\t",$0); print;}' |  awk '{ for (i=2;i<=NF;i+=2) $i="" } 1' | awk  '{gsub("  ","\t",$0); print;}' | awk  '{gsub(" ","",$0); print;}'   > $output/temp/${analname}_body_heading${outname}.Sfst
	awk '{gsub("=","\t",$0); print;}' $output/temp/${analname}_body${outname}.Sfst | awk '{ for (i=1;i<=NF;i+=2) $i="" } 1' | awk  '{gsub("  ","\t",$0); print;}' | awk  '{gsub(" ","",$0); print;}'  > $output/temp/${analname}_body2${outname}.Sfst
	cat $output/temp/${analname}_body_heading${outname}.Sfst $output/temp/${analname}_body2${outname}.Sfst > $output/temp/${analname}_prep${outname}.Sfst

	if  [[ -f ${output}/${analname}_avoidcols_${outname}.txt ]] ; then
		if [[ "$pops" != "NULL" ]] ; then
			COMBLIST2=$(awk '{print}' ORS=' ' ${output}/${analname}_avoidcols_${outname}.txt) 
			COMBLIST3=$(echo $COMBLIST2 | sed 's/[^ ][^ ]*/"&"/g')
				ZERO=$(echo $COMBLIST3 | sed -r 's/([^ ]+)/$i==\1/g') 
				TWO=$(echo $ZERO | sed -e 's/ /\  \|\|\ /g')
				ONE=$(echo " NR==1{for(i=1; i<=NF; i++) if (")
				THREE=$(echo ' ) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"} ')
				FOUR=$(echo ${ONE}${TWO}${THREE})
		else
			COMBLIST2=$(awk '{print}' ORS=' ' ${output}/${analname}_avoidcols_${outname}.txt) 
			COMBLIST3=$(echo $COMBLIST2 | sed 's/[^ ][^ ]*/"&"/g')
				ZERO=$(echo $COMBLIST3 | sed -r 's/([^ ]+)/$i!=\1/g') 
				TWO=$(echo $ZERO | sed -e 's/ /\  \&\&\ /g')
				ONE=$(echo " NR==1{for(i=1; i<=NF; i++) if (")
				THREE=$(echo ' ) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"} ')
				FOUR=$(echo ${ONE}${TWO}${THREE})
		fi
			#Use above expression to remove FST columns in similar populations
			awk "$FOUR" $output/temp/${analname}_prep${outname}.Sfst > $output/temp/${analname}_keep${outname}.Sfst
			NCOLZ=$(awk '{print NF}' $output/temp/${analname}_keep${outname}.Sfst | sort -nu | tail -n 1)
			if [ "$NCOLZ" -gt 1 ]; then
				awk '{if (NR!=1) {print}}' $output/temp/${analname}_keep${outname}.Sfst | awk '{ s = 0; for (i = 1; i <= NF; i++) s += $i; print $1, (NF > 1) ? s / (NF - 0) : 0; }' | awk '{print $2}' > $output/temp/${analname}_keep2${outname}.Sfst
			else
				awk '{if (NR!=1) {print}}' $output/temp/${analname}_keep${outname}.Sfst > $output/temp/${analname}_keep2${outname}.Sfst
			fi
			#Create final fst file for this analysis.
			paste $output/temp/${analname}_head${outname}.Sfst $output/temp/${analname}_keep2${outname}.Sfst | awk '{gsub("\t","",$0); print;}'  > $output/temp/${analname}_fstNA${outname}.Sfst
	else 
		declare -i NCOLZ=$(awk '{print NF}' $output/temp/${analname}_prep${outname}.Sfst | sort -nu | tail -n 1)
		echo $NCOLZ
		if [ "$NCOLZ" -gt 1 ]; then
			awk '{if (NR!=1) {print}}' $output/temp/${analname}_prep${outname}.Sfst | awk '{ s = 0; for (i = 1; i <= NF; i++) s += $i; print $1, (NF > 1) ? s / (NF - 0) : 0; }' | awk '{print $2}' > $output/temp/${analname}_keep2${outname}.Sfst
		else
			awk '{if (NR!=1) {print}}' $output/temp/${analname}_prep${outname}.Sfst > $output/temp/${analname}_keep2${outname}.Sfst
		fi
			#Create final fst file for this analysis.
			paste $output/temp/${analname}_head${outname}.Sfst $output/temp/${analname}_keep2${outname}.Sfst | awk '{gsub("\t","",$0); print;}'  > $output/temp/${analname}_fstNA${outname}.Sfst
	fi

	awk '!/na/' $output/temp/${analname}_fstNA${outname}.Sfst > $output/${analname}_${outname}.Sfst
	declare -i before=$(wc -l $output/${analname}_raw.Sfst |  cut -f1 -d' ')
	declare -i after=$(wc -l $output/${analname}_${outname}.Sfst | cut -f1 -d' ')
	declare -i loss=$(($before-$after))
	echo "ALERT: $loss SNPs were removed from FST analysis due to Ns, indels, or uninformative comparisons"
	echo "ALERT: $after SNPs remain"
	rm $output/temp/*${outname}*.Sfst*
fi

if [[ "$analtype" == "FET" ]]; then
	echo "ALERT: Subsetting FET file ${analname}_raw.fet"
	awk '$5 >= '$mincov' &&  $5  <= '$maxcov' && $4 >= '$propsnps' && $3 >= '$numsnps' ' $output/${analname}_raw.fet >  $output/temp/${analname}_reclassing${outname}.fet

				declare -i before=$(wc -l $output/${analname}_raw.fet |  cut -f1 -d' ')
				declare -i afterc=$(wc -l $output/temp/${analname}_reclassing${outname}.fet | cut -f1 -d' ')
				declare -i lossrc=$(($before-$afterc))
				echo "ALERT: $lossrc SNPs were removed from fet analysis due to new parameters"

	#Cut up fet file
	cut -d$'\t' -f 1-2 $output/temp/${analname}_reclassing${outname}.fet | gawk '$3=(FNR FS $3)' > $output/temp/${analname}_head${outname}.fet
	cut -d$'\t' -f 6- $output/temp/${analname}_reclassing${outname}.fet > $output/temp/${analname}_body${outname}.fet
	awk 'NR==1{print $0}' $output/temp/${analname}_body${outname}.fet | awk  '{gsub("=","\t",$0); print;}' |  awk '{ for (i=2;i<=NF;i+=2) $i="" } 1' | awk  '{gsub("  ","\t",$0); print;}' | awk  '{gsub(" ","",$0); print;}'   > $output/temp/${analname}_body_heading${outname}.fet
	awk '{gsub("=","\t",$0); print;}' $output/temp/${analname}_body${outname}.fet | awk '{ for (i=1;i<=NF;i+=2) $i="" } 1' | awk  '{gsub("  ","\t",$0); print;}' | awk  '{gsub(" ","",$0); print;}'  > $output/temp/${analname}_body2${outname}.fet
	cat $output/temp/${analname}_body_heading${outname}.fet $output/temp/${analname}_body2${outname}.fet > $output/temp/${analname}_prep${outname}.fet

	if  [[ -f ${output}/${analname}_avoidcols_${outname}.txt ]] ; then
		if [[ "$pops" != "NULL" ]] ; then
			COMBLIST2=$(awk '{print}' ORS=' ' ${output}/${analname}_avoidcols_${outname}.txt) 
			COMBLIST3=$(echo $COMBLIST2 | sed 's/[^ ][^ ]*/"&"/g')
				ZERO=$(echo $COMBLIST3 | sed -r 's/([^ ]+)/$i==\1/g') 
				TWO=$(echo $ZERO | sed -e 's/ /\  \|\|\ /g')
				ONE=$(echo " NR==1{for(i=1; i<=NF; i++) if (")
				THREE=$(echo ' ) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"} ')
				FOUR=$(echo ${ONE}${TWO}${THREE})
		else
			COMBLIST2=$(awk '{print}' ORS=' ' ${output}/${analname}_avoidcols_${outname}.txt) 
			COMBLIST3=$(echo $COMBLIST2 | sed 's/[^ ][^ ]*/"&"/g')
				ZERO=$(echo $COMBLIST3 | sed -r 's/([^ ]+)/$i!=\1/g') 
				TWO=$(echo $ZERO | sed -e 's/ /\  \&\&\ /g')
				ONE=$(echo " NR==1{for(i=1; i<=NF; i++) if (")
				THREE=$(echo ' ) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"} ')
				FOUR=$(echo ${ONE}${TWO}${THREE})
		fi
			#Use above expression to remove fet columns in similar populations
			awk "$FOUR" $output/temp/${analname}_prep${outname}.fet > $output/temp/${analname}_keep${outname}.fet
			awk '{if (NR!=1) {print}}' $output/temp/${analname}_keep${outname}.fet > $output/temp/${analname}_keep2${outname}.fet
			#Create final fet file for this analysis.
			
			rin=$output/temp/${analname}_keep2${outname}.fet
			rout=$output/temp/
			Rscript $BASEDIR/rscripts/r_combine_p.R $rin $rout
			

	else 
			awk '{if (NR!=1) {print}}' $output/temp/${analname}_keep${outname}.fet > $output/temp/${analname}_keep2${outname}.fet
			#Create final fet file for this analysis.
			
			rin=$output/temp/${analname}_keep2${outname}.fet
			rout=$output/temp/
			Rscript $BASEDIR/rscripts/r_combine_p.R $rin $rout
			
	fi
	
				paste $output/temp/${analname}_head${outname}.fet $output/temp/${analname}_keep2${outname}2.fet | awk '{gsub("\t","",$0); print;}'  > $output/temp/${analname}_fetNA${outname}.fet
				awk '!/na/' $output/temp/${analname}_fetNA${outname}.fet > $output/${analname}_${outname}.fet
				declare -i before=$(wc -l $output/${analname}_raw.fet |  cut -f1 -d' ')
				declare -i after=$(wc -l $output/${analname}_${outname}.fet | cut -f1 -d' ')
				declare -i loss=$(($before-$after))
				echo "ALERT: $after SNP positions remain"
				rm $output/temp/*${analname}*.fet*
fi

