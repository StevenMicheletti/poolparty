#!/bin/bash

#PoolParty_align_v0.76
#By Steven Micheletti

set -o pipefail

BASEDIR=$(dirname "$0")

if [[ ${BASEDIR} != *poolparty* ]];then
	SLOC=$(readlink -f "$0")
	BASEDIR=${SLOC%/*}
fi

echo $BASEDIR

if  ( [[ $(echo $1)  = "" ]] )  ; then
	echo "ERROR: You must provide a config file after the PPalign command"
	echo "for example, PPalign example_align.config"
	exit
fi
	echo "ALERT: $1 has been specified as the configuration file"

source $1

#Declare integers for numerical parameters
declare -i THREADZ
declare -i BQUAL
declare -i MINLENGTH
declare -i INWIN
declare -i MAPQ
declare -i SNPQ
declare -i MINDP

#Get date of run and save an additional config file with run parameters
RUNDATE=$(date +'%m_%d_%Y')

echo "ALERT: Beginning PoolParty run at $(date)"

########################################################################################################
####-----------------------------------------RUNNING-----------------------------------------------#####
########################################################################################################

	echo "ALERT: Performing error checking"

	#Check for samplelist file, quit if it does not exist
		if [[ ! -f $INDIR/"samplelist.txt" ]] ; then
 			 echo 'ERROR: File "samplelist.txt" is not there, aborting.'
   			 exit
		fi
	#Check samplelist for invalid characters
		if grep -q / $INDIR/"samplelist.txt" ; then
  			echo 'ERROR: "samplelist.txt" contains invalid characters. Do not include directory names or backslashes in samplelist'
			exit
		fi

	#Check for at least two entries in sample list (since pipeline input is paired-end data)
		leN=$(wc -l $INDIR/"samplelist.txt" | cut -f1 -d' ') 
		if  [[ $leN -lt 1 ]]  ; then
 			 echo 'ERROR: File "samplelist.txt" has less than two rows. Something wrong with file or not paired-end reads.'
   			 exit
		fi
	#Check that populations are designated in the samplelist
		PnUm=$(awk '{ a[$2]++ } END { for (b in a) { print b } }' $INDIR/"samplelist.txt")
		if [[ $PnUm  = "" ]] ; then
			echo 'ERROR: Something wrong with samplelist.txt . Are population numbers specified in column 2?'
			exit
		fi
		parray=($PnUm)
		parraynum=$(echo ${#parray[@]})
			echo "ALERT: You have specified $parraynum unique populations. If this sounds wrong then you done messed up..."
	#Check for genome fasta file
		if [[ ! -f ${GENOME} ]] ; then
 	 		echo "ERROR: ${GENOME} genome file is missing"
   		 	exit
		fi
	#Check for index file
		if [[ ! -f ${GENOME}.fai ]] ; then
 			 echo "ERROR:: ${GENOME}.fai index file is missing, place this in the same dir as the genome .fasta"
   			 exit
		fi
	#Check skip break is either on or off
		if [[ "$SPLITDISC" != "off" ]] && [[ "$SPLITDISC" != "on" ]] ; then
			echo "ERROR: Wrong Parameters: SPLITDISC must be either on or off"
			exit
		fi
	#Check individual analysis is either on or off
		if [[ "$INDCONT" != "off" ]] && [[ "$INDCONT" != "on" ]] ; then
			echo "ERROR: Wrong Parameters: INDCONT must be either on or off"
			exit
		fi
	#Check quality report is either on or off
		if [[ "$QUALREPORT" != "off" ]] && [[ "$QUALREPORT" != "on" ]] ; then
			echo "ERROR: Wrong Parameters: QUALREPORT must be either on or off"
			exit
		fi
	#Ensure base directories all exist
		if [ ! -d "$INDIR" ] || [ ! -d "$BBMAPDIR" ]  || [ ! -d "$POOL2" ] ; then
			echo "ERROR: Wrong Directories: Either INDIR, POOL2, or BBMAPDIR don't exist"
			exit
		fi
	#Ensure poolparty directory is intact
		if [ ! -d "$BASEDIR/rscripts" ] ; then
			echo "ERROR: Poolparty directory has been tampered with"
			exit
		fi
	#Ensure Popoolation2 directory structure has not been tampered with
		if [[ ! -f ${POOL2}/indel_filtering/"identify-indel-regions.pl" ]] ||  [[ ! -f ${POOL2}/indel_filtering/"filter-sync-by-gtf.pl" ]] \
			||  [[ ! -f ${POOL2}/"mpileup2sync.jar" ]]  ; then
			echo "ERROR: A component of Popoolation2 is missing. Ensure the directory structure in Popoolation2 is unaltered."
			echo "Ensure mpileup2sync.jar, filter-sync-by-gtf.pl, and identify-indel-regions.pl are all present"
			exit
		fi
	#Ensure BBMAP directory structure has not been tampered with
		if [[ ! -f ${BBMAPDIR}/"bbduk.sh" ]] ||  [[ ! -f ${BBMAPDIR}/resources/"adapters.fa" ]]  ; then
			echo "ERROR: A component of bbmap is missing. Ensure the directory structure in bbmap is unaltered."
			echo "Ensure that bbduk.sh and /resources/adapters.fa exist in the bbmap directory"
			exit
		fi
	#Ensure that all dependencies exist on the system (both the path and the command have to fail)
		if  ( [[ $(command -v "$FASTQC")  = "" ]]  &&  [[ ! -f $FASTQC ]] ) || ( [[ $(command -v "$BWA")  = "" ]]  &&  [[ ! -f $BWA ]] ) \
				|| ( [[ $(command -v "$SAMTOOLS")  = "" ]]  &&  [[ ! -f $SAMTOOLS ]] ) ||  ( [[ $(command -v "$PICARDTOOLS")  = "" ]]  &&  [[ ! -f $PICARDTOOLS ]] ) \
				|| ( [[ $(command -v "$SAMBLASTER")  = "" ]]  &&  [[ ! -f $SAMBLASTER ]] ) || ( [[ $(command -v "$BCFTOOLS")  = "" ]]  &&  [[ ! -f $BCFTOOLS ]] ) ; then
			echo "ERROR: One or more dependencies are incorrect, double check dependency locations and names"
			echo "Type each dependency into the terminal, as it is listed in the config file, and ensure that it initiates"		
			exit
		fi
	#Make sure java memory is set properly
		if ! printf '%s\n' "$KMEM" | grep -Fqe "Xmx"; then
			echo "ERROR: Java memory set incorrectly, check value and ensure it followes Xmx#g format"
			exit
		fi
	#Ensure other parameters are within a valid range
		if  [[ ! "$THREADZ" =~ ^[0-9]+$ ]] || [[ ! "$BQUAL" =~ ^[0-9]+$ ]] || [[ ! "$MINLENGTH" =~ ^[0-9]+$ ]] || [[ ! "$INWIN" =~ ^[0-9]+$ ]] \
				|| [[ ! "$MAPQ" =~ ^[0-9]+$ ]] || [[ ! "$MINDP" =~ ^[0-9]+$ ]] ;  then
        		echo "ERROR: THREADZ, BQUAL, MINLENGTH, INWIN, and MAPQ must all be positive integers. Check parameters"
			exit
		fi
		if  [[ ! "$MAF" =~ ^[+-]?[0-9]+\.?[0-9]*$  ]] ; then 
				echo "ERROR: MAF must be a value between 0 and 1"
				exit
		fi
	#Check for R, perl, java on the system
		if  [[ $(command -v java)  = "" ]] ; then
			echo "ERROR: java not detected on system. java is required for multiple packages"
			echo "java should initiate when 'java' is typed in terminal"
			exit
		fi
		if  [[ $(command -v perl)  = "" ]] ; then
			echo "ERROR: perl not detected on system. perl is required for Popoolation2"
			echo "perl should initiate when 'perl' is typed in terminal"
			exit
		fi
		if  [[ $(command -v Rscript)  = "" ]] ; then
			echo "ERROR: R not detected on system. R is required for analysis "
			echo "R should initiate when 'Rscript' is typed in terminal"
			exit
		fi
	#Check for piping ability 
		if  [[ $(command -v mkfifo)  = "" ]] ; then
			echo "WARNING: piping 'mkfifo' not detected on system or available on drive. This may cause issues in downstream analyses"
			echo "Edit PPalign.sh and redirect mkfifo to another drive"
		fi
	#check for gawk
		if  [[ $(command -v gawk)  = "" ]] ; then
			echo "ERROR: gawk not detected on system. gawk is usually standard on Linux systems. Install then retry"
			exit
		fi
	#Check for process substitution 
		if  [[ $(command -v cat <(date); echo $? )  = "" ]] ; then
			echo "WARNING: process substitution does not appear to be working on this system. Errors may arise"
		fi
	#Load R to check for dependencies 
		if  [[ -f $OUTDIR/R_ERROR.txt ]] ; then
			rm $OUTDIR/R_ERROR.txt 
		fi		

		Rscript $BASEDIR/rscripts/r_align_check.R $OUTDIR

		if  [[ -f $OUTDIR/R_ERROR.txt ]] ; then
			echo "ERROR: R dependency check failed, install dependencies manually"
			exit
		fi

	echo "ALERT: Parameter check passed. moving on..."

#########################################################################################################
#Make the output directory if it does not exist 
	if [ ! -d "$OUTDIR/" ]; then
		mkdir $OUTDIR/
	fi

#Copy configuration file so you know what you did later
	cp $1 ${OUTDIR}/${OUTPOP}_${RUNDATE}.config

#Remove any potential extra characters (/r) and sort  samplelist.txt
	awk '{ sub("\r$", ""); print $0 }'  ${INDIR}/"samplelist.txt" > ${INDIR}/samplelist.tmp && mv ${INDIR}/samplelist.tmp ${INDIR}/samplelist.txt
	sort -f ${INDIR}/"samplelist.txt" | awk '{print $1}' | grep -P -v '^\s*$' >  $OUTDIR/${OUTPOP}_sample_files.txt
	sort -f ${INDIR}/"samplelist.txt" | awk '{print $2}' | grep -P -v '^\s*$' > $OUTDIR/${OUTPOP}_sample_pops.txt
	FILES=$(cat $OUTDIR/${OUTPOP}_sample_files.txt)
	POPS=$(cat $OUTDIR/${OUTPOP}_sample_pops.txt)
		aCHCK1=$(wc -l $OUTDIR/${OUTPOP}_sample_files.txt | awk '{print $1}' )
		aCHCK2=$(wc -l $OUTDIR/${OUTPOP}_sample_pops.txt | awk '{print $1}' )
			#Ensure each sample has been given a population designation
			if [[ $aCHCK1 != $aCHCK2  ]] ; then
				echo "ERROR: Something wrong with samplelist.txt. $aCHCK1 samples but only $aCHCK2 population designations.."
				echo "Check for empty lines in samplelist.txt"
				exit
			fi

#Check mate pairs; fail if order is incorrect
	cut -d_ -f1  $OUTDIR/${OUTPOP}_sample_files.txt >  $OUTDIR/${OUTPOP}_prefixes.txt

	CHECK=$(cat $OUTDIR/${OUTPOP}_prefixes.txt)
	CHECK2=($CHECK)

	for (( i=0; i<${#CHECK2[@]} ; i+=2 )) ; do
		echo  Checking pairs "${CHECK2[i+1]}" and "${CHECK2[i]}"
		if [ "${CHECK2[i+1]}" == "${CHECK2[i]}" ]; then
		echo "ALERT: Correct Mates, proceeding"
		single=off
	else
		echo "ALERT: Naming convention suggests single-end reads. Proceeding using single-end fqs."
		echo "If reads are paired-end, check naming convention and try again."
		single=on
		break
	fi
	done

if [[ "$single" =~(off)$ ]] ; then
#Create prefix file; this indicates prefix before the first "_" which should be library identifiers
	cut -d_ -f1 $OUTDIR/${OUTPOP}_sample_files.txt | awk '!seen[$0]++' | sed '/^$/d'  >  $OUTDIR/${OUTPOP}_prefixes.txt
	printf "ALERT: Using libraries:\n$(cat $OUTDIR/${OUTPOP}_prefixes.txt)\n"
		declare -i CHCK1=$(wc -l $OUTDIR/${OUTPOP}_prefixes.txt | cut -f1 -d' ' )
		declare -i CHCK2=$(wc -l $OUTDIR/${OUTPOP}_sample_files.txt | cut -f1 -d' '  )
		declare -i CHCK3=($CHCK2)/2
		if [[ $CHCK1 != $CHCK3  ]] ; then
        	echo "ERROR: number of samples and prefixes don't match up ($CHCK1 vs $CHCK3) ; there is likely something wrong with the filenames or samplelist"
			echo "Note that file name prefixes (text before the first underscore) must be unique for the paired files!"
			echo "Check _sample_files.txt, paired-end libraries should be stacked on top of one another; if this isn't the case change the naming convention"
			exit
		fi
fi

if [[ "$single" =~(on)$ ]] ; then
#Create prefix file; this indicates prefix before the first "_" which should be library identifiers
	cut -d_ -f1 $OUTDIR/${OUTPOP}_sample_files.txt | awk '!seen[$0]++' | sed '/^$/d'  >  $OUTDIR/${OUTPOP}_prefixes.txt
	printf "ALERT: Using single-ended libraries:\n$(cat $OUTDIR/${OUTPOP}_prefixes.txt)\n"
fi

#Get number of populations that are represented by the libraries in your samplelist
	if  [[ -f $OUTDIR/${OUTPOP}_poplist.txt ]] ; then
		rm  $OUTDIR/${OUTPOP}_poplist.txt 
	fi

	array=($FILES)
	array2=($POPS)
	echo "ALERT: Checking ${#array[@]} files in samplelist" 
	if [[ "$single" =~(off)$ ]] ; then
		for (( i=0; i<${#array[@]} ; i+=2 )) ; do
			e=${array[i]%%.*}
			f=${array2[i]}
			oZ=$(echo $e $f)
			echo ${oZ} >> $OUTDIR/${OUTPOP}_poplist.txt
			#remove any introduced characters
			awk '{ sub("\r$", ""); print $0 }'  $OUTDIR/${OUTPOP}_poplist.txt > $OUTDIR/${OUTPOP}_poplist.tmp && mv $OUTDIR/${OUTPOP}_poplist.tmp $OUTDIR/${OUTPOP}_poplist.txt
		done
	fi
	if [[ "$single" =~(on)$ ]] ; then
		for (( i=0; i<${#array[@]} ; i+=1 )) ; do
			e=${array[i]%%.*}
			f=${array2[i]}
			oZ=$(echo $e $f)
			echo ${oZ} >> $OUTDIR/${OUTPOP}_poplist.txt
			#remove any introduced characters
			awk '{ sub("\r$", ""); print $0 }'  $OUTDIR/${OUTPOP}_poplist.txt > $OUTDIR/${OUTPOP}_poplist.tmp && mv $OUTDIR/${OUTPOP}_poplist.tmp $OUTDIR/${OUTPOP}_poplist.txt
		done
	fi
		rm $OUTDIR/${OUTPOP}_sample_pops.txt
		rm $OUTDIR/${OUTPOP}_sample_files.txt
		rm $OUTDIR/${OUTPOP}_prefixes.txt 

#Get anchored chromosome lengths in bp. Prints a copy to the rundir for later analyses
	echo "ALERT: Getting genome anchored index..."
		HNAM=$(cut -f1 ${GENOME}.fai |  sed -re "s/[^a-zA-Z]*//g"  |  awk '!x[$0]++')
		HNAM2=$(echo $HNAM)
		echo "ALERT: $HNAM2 are your chromosome and/or scaffold headings"
		if [ -z "$SCAHEAD" ]; then
			cut -f1,2 ${GENOME}.fai >  $OUTDIR/${OUTPOP}_CHRbp1.txt
		else
			cut -f1,2 ${GENOME}.fai  | grep -v "^${SCAHEAD}"  >  $OUTDIR/${OUTPOP}_CHRbp1.txt
		fi
	awk '{print $1}' $OUTDIR/${OUTPOP}_CHRbp1.txt  | awk  '$2="1"' | awk '{gsub(" ","\t",$0); print;}' > $OUTDIR/${OUTPOP}_CHRbp2.txt
	cat $OUTDIR/${OUTPOP}_CHRbp1.txt $OUTDIR/${OUTPOP}_CHRbp2.txt > $OUTDIR/${OUTPOP}_CHRbp.txt
	rm $OUTDIR/${OUTPOP}_CHRbp1.txt ; rm $OUTDIR/${OUTPOP}_CHRbp2.txt

		if  [[ -f $OUTDIR/${OUTPOP}_names.txt ]] ; then
			rm  $OUTDIR/${OUTPOP}_names.txt 
		fi

# Make separate files for pops, indicating which library belongs to which population
	if [ ! -d "$OUTDIR/pops" ]; then
		mkdir $OUTDIR/pops
	fi

#Split into one file per population to combine
	awk '{print $1 > "'${OUTDIR}/pops/pop_'"$2".txt"}' $OUTDIR/${OUTPOP}_poplist.txt
	popfiles=$(awk '{ a[$2]++ } END { for (b in a) { print b } }' $INDIR/"samplelist.txt")
	printf '%s\n' "${popfiles[@]}" | awk '{print "pop_" $0;}' | sort > ${OUTDIR}/pops/${OUTPOP}_files_for_pops.txt
	awk '{ sub("\r$", ""); print $0 }'  ${OUTDIR}/pops/${OUTPOP}_files_for_pops.txt > ${OUTDIR}/pops/${OUTPOP}_files_for_pops.tmp && mv ${OUTDIR}/pops/${OUTPOP}_files_for_pops.tmp ${OUTDIR}/pops/${OUTPOP}_files_for_pops.txt


#Trim by quality score, uses BBMAP (java-based)
##! ADDITIONAL TRIM PARAMETERS CAN BE ADDED BELOW AFTER "${BBMAPDIR}/bbduk.sh" !##
##CHECK BBMAP DOCUMENTATION FOR ADDITIONAL OPTIONS AND MODIFY LINE BELOW ##
	if [ ! -d "$OUTDIR/trimmed" ]; then
		mkdir $OUTDIR/trimmed
	fi

	array=($FILES)
	echo trimming  ${#array[@]} files 

	if [[ "$single" =~(off)$ ]] ; then
		for (( i=0; i<${#array[@]} ; i+=2 )) ; do
			b=${array[i]%%.*}
			c=${array[i+1]%%.*}
			o1=${OUTDIR}/trimmed/${b}.trim_1
			o2=${OUTDIR}/trimmed/${b}.trim_2 

			if  [[ -f $o1 ]] && [[ -f $o2 ]]  ; then
				echo "ALERT: $o1 pair exists; skipping"
				echo ${b} >> $OUTDIR/${OUTPOP}_names.txt
			fi

			if [[ ! -f $o1 ]] && [[ ! -f $o2 ]] ; then
			#Add modifications below if needed#
				echo ${b} >> $OUTDIR/${OUTPOP}_names.txt
				if [[ -f ${OUTDIR}/BAM/${b}_filtered.bam ]] ; then
					echo "ALERT: ${b}_filtered.bam exists and will not be re-trimmed or re-aligned."
				else
					${BBMAPDIR}/bbduk.sh -${KMEM} in=$INDIR/${array[i]} in2=$INDIR/${array[i+1]} out=$o1 out2=$o2 \
						ref=${BBMAPDIR}/resources/adapters.fa ktrim=r k=23 kmask=f tpe tbo qtrim=r trimq=${BQUAL} minlength=${MINLENGTH} minavgquality=${BQUAL} threads=${THREADZ} stats=$OUTDIR/trimmed/${b}_trimstats.txt
				fi
			fi
		done
	fi
	
	if [[ "$single" =~(on)$ ]] ; then
		for (( i=0; i<${#array[@]} ; i+=1 )) ; do
			b=${array[i]%%.*}
			c=${array[i+1]%%.*}
			o1=${OUTDIR}/trimmed/${b}.trim_1

			if  [[ -f $o1 ]] ; then
				echo "ALERT: $o1  exists; skipping"
				echo ${b} >> $OUTDIR/${OUTPOP}_names.txt
			fi

			if [[ ! -f $o1 ]] ; then
				echo ${b} >> $OUTDIR/${OUTPOP}_names.txt
				if [[ -f ${OUTDIR}/BAM/${b}_filtered.bam ]] ; then
					echo "ALERT: ${b}_filtered.bam exists and will not be re-trimmed or re-aligned."
				else
					#Add modifications below if needed#
					${BBMAPDIR}/bbduk.sh -${KMEM} in=$INDIR/${array[i]} out=$o1  \
					ref=${BBMAPDIR}/resources/adapters.fa ktrim=r k=23 kmask=f tpe tbo qtrim=r trimq=${BQUAL} minlength=${MINLENGTH} minavgquality=${BQUAL} threads=${THREADZ} stats=$OUTDIR/trimmed/${b}_trimstats.txt
				fi
			fi
		done
	fi

		#remove weird characters from names
		awk '{ sub("\r$", ""); print $0 }'  $OUTDIR/${OUTPOP}_names.txt > $OUTDIR/${OUTPOP}_names.tmp && mv $OUTDIR/${OUTPOP}_names.tmp $OUTDIR/${OUTPOP}_names.txt
	
	#Check that trimmed folder was written to
	if [ -z "$(ls -A ${OUTDIR}/trimmed)" ]; then
		echo "WARNING: No trimmed files were produced"
	fi

	if [[ "$QUALREPORT" =~(on)$ ]] && [ ! -z "$(ls -A ${OUTDIR}/trimmed)" ] ; then
	# Quality report of trimmed files (FASTQC). This runs in the background as alignments are produced
		if [ ! -d "$OUTDIR/quality" ]; then
			mkdir $OUTDIR/quality
		fi
			echo "ALERT: FASTQC has started running in background"
			for i in ${OUTDIR}/trimmed/*.{trim_1,trim_2}; do
				ni=$(echo $i | awk -F/ '{print $NF}') 
				if [[ ! -f ${OUTDIR}/quality/${ni}_fastqc.zip ]] ; then
					${FASTQC} -q -o ${OUTDIR}/quality $i 
				fi
			done &
	fi

#Alignment and duplicate removal (BWA, SAMBLASTER, SAMTOOLS)
##!ADDITIONAL ALIGNMENT PARAMETERS CAN BE ADDED BELOW AFTER "${BWA} mem" !##
#CHECK bwa mem DOCUMENTATION FOR ADDITIONAL OPTIONS AND MODIFY LINE BELOW#
	
	#Make folder for BAMS
	echo "ALERT: Beginning BWA mem at $(date)"
	
	ITER=1
	if [ ! -d "$OUTDIR/BAM" ]; then
	mkdir $OUTDIR/BAM
	fi
	
	for i in $(cat $OUTDIR/${OUTPOP}_names.txt); do
		if [[ -f ${OUTDIR}/BAM/${i}_aligned.bam ]] || [[ -f ${OUTDIR}/BAM/${i}_sorted.bam ]] || [[ -f ${OUTDIR}/BAM/${i}_filtered.bam ]]; then
			echo "ALERT: ${i}_aligned.bam or ${i}_sorted.bam exists or ${i}_filtered.bam; skipping"
		else
			if [[ "$single" =~(off)$ ]] ; then
				if [[ "$SPLITDISC" =~(on)$ ]] ; then 
					echo "ALERT: Aligning with discordant and split read production"
					#Makes unique lane and pop name for each library
					a="@RG\tID:SMP"
					b="\tPL:ILLUMINA\tLB:Pooled\tSM:"
					c=$a$ITER$b$i
					nice -n 19 ${BWA} mem -M -t $THREADZ -R $c \
					"${GENOME}" \
					"${OUTDIR}/trimmed/$i.trim_1" "${OUTDIR}/trimmed/$i.trim_2" | ${SAMBLASTER} -M -r -d "${OUTDIR}/BAM/$i.disc.sam" -s "${OUTDIR}/BAM/$i.split.sam"  | ${SAMTOOLS} view -Sb -q ${MAPQ} - > "${OUTDIR}/BAM/${i}_aligned.bam"
				fi	
			
				if [[ "$SPLITDISC" =~(off)$ ]] ; then 
					echo "ALERT: Aligning without discordant and split read production"
					#Makes unique lane and pop name for each library
					a="@RG\tID:SMP"
					b="\tPL:ILLUMINA\tLB:Pooled\tSM:"
					c=$a$ITER$b$i
					nice -n 19 ${BWA} mem -M -t $THREADZ -R $c \
					"${GENOME}" \
					"${OUTDIR}/trimmed/$i.trim_1" "${OUTDIR}/trimmed/$i.trim_2" | ${SAMBLASTER} -M -r | ${SAMTOOLS} view -Sb -q ${MAPQ} - > "${OUTDIR}/BAM/${i}_aligned.bam"
				fi
			fi
			if [[ "$single" =~(on)$ ]] ; then
				if [[ "$SPLITDISC" =~(on)$ ]] ; then 
					echo "ALERT: Cannot perform discordant/split-end analyses on single-end reads! ignoring"
				fi	
				echo "ALERT: Aligning without discordant and split read production"
				#Makes unique lane and pop name for each library
				a="@RG\tID:SMP"
				b="\tPL:ILLUMINA\tLB:Pooled\tSM:"
				c=$a$ITER$b$i
				nice -n 19 ${BWA} mem -M -t $THREADZ -R $c \
				"${GENOME}" \
				"${OUTDIR}/trimmed/$i.trim_1" | ${SAMBLASTER} -M -r | ${SAMTOOLS} view -Sb -q ${MAPQ} - > "${OUTDIR}/BAM/${i}_aligned.bam"
			fi
		fi

	ITER="$(($ITER + 1))"
	done

	echo "ALERT: Finished BWA mem at $(date)"
		#reports will have read alignment results for each bam file
		if [ ! -d "$OUTDIR/reports" ]; then
			mkdir $OUTDIR/reports
		fi

# Sorting BAM files and filtering unpaired reads (PICARDTOOLS, SAMTOOLS)
	for i in $(cat $OUTDIR/${OUTPOP}_names.txt); do
		if  [[ -f ${OUTDIR}/BAM/${i}_sorted.bam ]] || [[ -f ${OUTDIR}/BAM/${i}_filtered.bam ]] ; then
			echo "ALERT: ${i}_sorted.bam or filtered.bam exists; skipping"
		else
			echo "ALERT: Picardtools started sorting ${OUTDIR}/BAM/${i}_aligned.bam at $(date) "
			nice -n 19 java -XX:ParallelGCThreads=$THREADZ -${KMEM} -Djava.io.tmpdir=`pwd`/tmp -jar ${PICARDTOOLS} SortSam I= ${OUTDIR}/BAM/${i}_aligned.bam O= ${OUTDIR}/BAM/${i}_sorted.bam VALIDATION_STRINGENCY=SILENT QUIET=true SO=coordinate TMP_DIR=`pwd`/tmp 
		fi
	done
		echo "ALERT: Picardtools finished sorting all BAMS at $(date) " 
	#Making reports for sorted BAM files
		echo "ALERT: Alignment reports started at $(date) and will run in the background"
		for i in $(cat $OUTDIR/${OUTPOP}_names.txt); do
			if [[ ! -f ${OUTDIR}/reports/${i}_${OUTPOP}_aln_report.txt ]] ; then
				nice -n 19 ${SAMTOOLS} flagstat ${OUTDIR}/BAM/${i}_sorted.bam > ${OUTDIR}/reports/${i}_${OUTPOP}_aln_report.txt
			fi
		done &

#Create semaphore to parallel run at specified 'THREADZ'
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

#Filtering bams; removing junk. Runs in parallel.
##! ADDITIONAL BAM FILTER PARAMETERS CAN BE ADDED BELOW AFTER "${SAMTOOLS} view " !##
		task(){
			if [[ -f ${OUTDIR}/BAM/${i}_filtered.bam ]] ; then
				echo "ALERT: ${i}_filtered.bam exists; skipping"
			else
				if [[ "$single" =~(off)$ ]] ; then
					echo "ALERT: samtools is filtering ${OUTDIR}/BAM/${i}_sorted.bam at $(date) "
					nice -n 19 ${SAMTOOLS} view -b -F 0x04 -f 0x02 ${OUTDIR}/BAM/${i}_sorted.bam > ${OUTDIR}/BAM/${i}_filtered.bam 
				else
					cp ${OUTDIR}/BAM/${i}_sorted.bam ${OUTDIR}/BAM/${i}_filtered.bam 
				fi
			fi
			if [[ ! -f ${OUTDIR}/BAM/${i}_filtered.bam ]] && [[ "$single" =~(on)$ ]] ; then

				echo "ALERT: samtools is filtering ${OUTDIR}/BAM/${i}_sorted.bam at $(date) "
				nice -n 19 ${SAMTOOLS} view -b ${OUTDIR}/BAM/${i}_sorted.bam > ${OUTDIR}/BAM/${i}_filtered.bam 
			fi
			}

		N=$THREADZ
		open_sem $N
			for i in $(cat $OUTDIR/${OUTPOP}_names.txt); do
				run_with_lock task $i
			done ; wait

		echo "ALERT: Samtools finished filtering all BAMs at $(date) "

	#Remove initial bams to reduce storage (pretty useless at this point)
	for i in $(cat $OUTDIR/${OUTPOP}_names.txt); do
		if [[ -f ${OUTDIR}/BAM/${i}_aligned.bam ]] ; then
			rm ${OUTDIR}/BAM/${i}_aligned.bam
		fi
	done	


	if [ ! -z $ALIGNONLY ] &&  [[ "$ALIGNONLY" =~(on)$ ]] ; then
			echo "ALERT: PPalign Alignment Only completed at $(date) "
			exit
	fi

#Combine BAMS into specified populations. Uses parallel
			POPZ=$(cat ${OUTDIR}/pops/${OUTPOP}_files_for_pops.txt)	
			declare -a farray=($POPZ)
			echo "ALERT: Samtools is combining ${#farray[@]} populations into BAMS at $(date)"

task() {
			declare -i pnuM=$(wc -l ${OUTDIR}/pops/${i}.txt |  cut -f1 -d' ')
			declare -a barray=$(awk '{print "'${OUTDIR}/BAM/'"$0"_filtered.bam"}' ${OUTDIR}/pops/${i}.txt)
			
			if [[ -f ${OUTDIR}/BAM/${i}.bam ]] ; then
				echo "ALERT: Population ${i}.bam exists; skipping"
			else

				#If population consists of one library, simply duplicate the bam file. 
				if [ $pnuM -lt 2 ] ; then 
					cp $barray ${OUTDIR}/BAM/${i}.bam
				else
					echo "ALERT: Samtools is merging $pnuM bams into populations ${i}.bam at $(date) "
					samtools merge -f ${OUTDIR}/BAM/${i}.bam $barray
				fi
			fi
	}

	N=$THREADZ
	open_sem $N
		for i in "${farray[@]}" ; do
			run_with_lock task $i
	done ; wait
	echo "ALERT: Samtools finished merging all bams at $(date) "

#Call SNPs and print variant sites
	echo "ALERT: bcftools is calling variants  at $(date) this might take awhile..."
	COMBINE=$(sed 's,^,'$OUTDIR/BAM/', ; s/$/.bam/' $OUTDIR/pops/${OUTPOP}_files_for_pops.txt)
	echo "ALERT: See ${OUTPOP}_files_for_pops.txt for population merge order"
	
	 if [[ -f $OUTDIR/${OUTPOP}_variants.txt ]] ; then
		echo "ALERT: ${OUTPOP}_variants.txt already exists, skipping variant calling "
	 else

	#Make the filter directory if it does not exist 
		if [ ! -d "$OUTDIR/filters" ]; then
			mkdir $OUTDIR/filters
		fi
		
		if [ ! -z $MAKEVCF ]; then
			${BCFTOOLS} view  -i  'MIN(DP)>'$MINDP' & MIN(QUAL)>'$SNPQ' ' $MAKEVCF  > $OUTDIR/${OUTPOP}_Qualtemp.VCF
				declare -i after=$(grep -v "^#" $OUTDIR/${OUTPOP}_Qualtemp.VCF | wc -l) 
				declare -i lost=$(($before-$after))
				echo "ALERT: $lost SNPs removed due to QUAL < $SNPQ and total DP < $MINDP "
			${BCFTOOLS} view  -i  'MAF[0]> '$MAF' ' $OUTDIR/${OUTPOP}_Qualtemp.VCF > $OUTDIR/${OUTPOP}.VCF
				declare -i afterm=$( grep -v "^#" $OUTDIR/${OUTPOP}.VCF | wc -l) 
				declare -i lostm=$(($after-$afterm))
				echo "ALERT: Additional $lostm SNPs removed due to global MAF < $MAF "
				echo "ALERT: $afterm total SNPs retained after SNP calling"
				grep -v "^#"  $OUTDIR/${OUTPOP}.VCF | awk '{gsub("=|;", "\t", $0); print;}' \
						| awk '{print $1,$2,$6,$9}'  >  $OUTDIR/${OUTPOP}_variants.txt
		else
			${SAMTOOLS} mpileup -uf ${GENOME} -B $COMBINE |  ${BCFTOOLS} call --threads ${THREADZ} -mv -Ov >  $OUTDIR/${OUTPOP}_full.VCF 
				declare -i before=$(grep -v "^#" $OUTDIR/${OUTPOP}_full.VCF | wc -l)
				echo "ALERT: $before SNPs total SNPS called without filters"
			${BCFTOOLS} view  -i  'MIN(DP)>'$MINDP' & MIN(QUAL)>'$SNPQ' ' $OUTDIR/${OUTPOP}_full.VCF  > $OUTDIR/${OUTPOP}_Qualtemp.VCF
				declare -i after=$(grep -v "^#" $OUTDIR/${OUTPOP}_Qualtemp.VCF | wc -l) 
				declare -i lost=$(($before-$after))
				echo "ALERT: $lost SNPs removed due to QUAL < $SNPQ and total DP < $MINDP "
			${BCFTOOLS} view  -i  'MAF[0]> '$MAF' ' $OUTDIR/${OUTPOP}_Qualtemp.VCF > $OUTDIR/${OUTPOP}.VCF
				declare -i afterm=$( grep -v "^#" $OUTDIR/${OUTPOP}.VCF | wc -l) 
				declare -i lostm=$(($after-$afterm))
				echo "ALERT: Additional $lostm SNPs removed due to global MAF < $MAF "
				echo "ALERT: $afterm total SNPs retained after SNP calling"
				grep -v "^#"  $OUTDIR/${OUTPOP}.VCF | awk '{gsub("=|;", "\t", $0); print;}' \
						| awk '{print $1,$2,$6,$9}'  >  $OUTDIR/${OUTPOP}_variants.txt
		fi

		if [[ ! -s $OUTDIR/${OUTPOP}_full.VCF ]] ; then
			echo "ERROR: Variants were not filtered, $OUTDIR/${OUTPOP}_variants.txt was not produced. Check for sufficient memory"
			exit
		fi
			rm $OUTDIR/${OUTPOP}_Qualtemp.VCF
			rm $OUTDIR/${OUTPOP}_full.VCF
		
		awk '!/IDV|,/' $OUTDIR/${OUTPOP}_variants.txt > $OUTDIR/filters/${OUTPOP}_SNP_variants.txt
		awk '/IDV/'  $OUTDIR/${OUTPOP}_variants.txt > $OUTDIR/filters/${OUTPOP}_indel_variants.txt
				numSNP=$(wc -l $OUTDIR/filters/${OUTPOP}_SNP_variants.txt|  cut -f1 -d' ')
				numINDEL=$(wc -l $OUTDIR/filters/${OUTPOP}_indel_variants.txt |  cut -f1 -d' ')

		echo "ALERT: Of the remaining SNPs, there are $numSNP SNPs and $numINDEL INDels"
		echo "ALERT: Variant calling and filtering done at $(date) "
	fi
	
	if [[ -f  $OUTDIR/${OUTPOP}.mpileup ]] ; then
		echo "ALERT: ${OUTPOP}.mpileup already exists, skipping mpileup creation calling "
	 else

		echo "ALERT: mpileups are being created at  $(date)"
		
			#Get Columns of mpileup and determine population number
			declare -i NCOL=$(($parraynum * 3 + 3))
			# Start at 3, increment by NPOPS until NCOL, i.e, get column of each pop
			POPSEQ=$(seq 4 3 $NCOL)
			#Select columns to run analyses on 
			declare -a arr2=($POPSEQ)
			#Subset
			ONE=$(echo $POPSEQ | sed -r 's/([^ ]+)/$\1/g')
			TWO=$(echo $ONE | awk '$NF=$NF "}" ' OFS=",")
			THREE='{print $1,$2,$3,'
			FOUR=$(echo ${THREE}${TWO})	
			${SAMTOOLS} mpileup -f ${GENOME} -B $COMBINE | awk "$FOUR"  > $OUTDIR/${OUTPOP}_stats.mpileup &  PIDMIX=$!
			gawk 'NR==FNR{a[$1,$2]=$5;next} ($1,$2) in a{print $0, a[$1,$2]}' <(awk '{print $0}' $OUTDIR/${OUTPOP}_variants.txt) <(${SAMTOOLS} mpileup -f ${GENOME} -B $COMBINE)| awk '{gsub(" ","",$0); print;}' > $OUTDIR/${OUTPOP}.mpileup &  PIDIOS=$!
			wait $PIDIOS
			wait $PIDMIX
			echo "ALERT: Mpileups created at $(date)"
	fi
	
		if [[ ! -s $OUTDIR/${OUTPOP}.mpileup ]]  ; then
			echo "ERROR: Mpileup is empty. Check for sufficient memory during samtools mpileup"
			exit
		fi
		
	#Filter by indels and create sync format
		echo "ALERT: Identifying indel regions and creating sync format at $(date)"
			if [[ -f $OUTDIR/${OUTPOP}.sync ]] ; then
				echo "ALERT: Filtered sync file already exists, skipping this step "
			else
				perl ${POOL2}/indel_filtering/identify-indel-regions.pl --input $OUTDIR/${OUTPOP}.mpileup --indel-window $INWIN --output ${OUTDIR}/${OUTPOP}.gtffile --min-count 1
				nice -n 19 java -ea -${KMEM} -Djava.io.tmpdir=`pwd`/tmp -jar ${POOL2}/mpileup2sync.jar --min-qual 1 --fastq-type sanger --input $OUTDIR/${OUTPOP}.mpileup --output $OUTDIR/${OUTPOP}_temp.sync --threads $THREADZ 
				perl ${POOL2}/indel_filtering/filter-sync-by-gtf.pl --input ${OUTDIR}/${OUTPOP}_temp.sync --gtf ${OUTDIR}/${OUTPOP}.gtffile --output ${OUTDIR}/${OUTPOP}.sync
			echo "ALERT: Done identifying indel regions and creating sync format at $(date)"
				declare -i before=$(wc -l < $OUTDIR/${OUTPOP}_temp.sync )
				declare -i after=$(wc -l < $OUTDIR/${OUTPOP}.sync)
				declare -i lost="$(($before - $after))"
				lostp=$((100-100*$after/$before))
				echo "ALERT: With an indel window of $INWIN bp you lost $lost SNPs or $lostp % " 
				rm $OUTDIR/${OUTPOP}_temp.sync 
			fi				
			
		if [[ ! -s $OUTDIR/${OUTPOP}.sync  ]]  ; then
			echo "ERROR: Sync is empty. Check for sufficient memory during popoolation2 pl scripts"
			exit
		fi
	
# For each population listed in _files_for_pops, combined the filtered bam for individuals in that file as mpileups, while only keeping variants With parrelization
if [[ "$INDCONT" =~(on)$ ]] ; then

		if [ ! -d "$OUTDIR/inds" ]; then
			mkdir $OUTDIR/inds
		fi
		echo "ALERT: Creating individual mpileups contribution in ${#farray[@]} populations at $(date)"

	task() {
		if [[ -f $OUTDIR/inds/${file}.mpileup ]] ; then
			echo "ALERT: ${file}.mpileup already exists, skipping individual mpileup "
		else
			#Combine all individuals from same population into same mpileup
			COMBINE=$(sed 's,^,'$OUTDIR/BAM/', ; s/$/_filtered.bam/' ${OUTDIR}/pops/$file.txt)
			gawk 'NR==FNR{a[$1,$2]=$5;next} ($1,$2) in a{print $0, a[$1,$2]}' <(awk '{print $0}' ${OUTDIR}/${OUTPOP}_variants.txt) <(${SAMTOOLS} mpileup -aa -f ${GENOME} -B $COMBINE) | awk '{gsub(" ","",$0); print;}' > $OUTDIR/inds/${file}.mpileup
		fi
	}

		N=$THREADZ
		open_sem $N
			for file in $(cat $OUTDIR/pops/${OUTPOP}_files_for_pops.txt) ; do
				run_with_lock task $i
		done ; wait
		echo "ALERT: Individual mpileups created at $(date) "

	task() {
			if [[ -f $OUTDIR/inds/${file}_RindSTATs.rin ]] ; then
				echo "ALERT: ${file}_RindSTATs.rin exists, skipping this file "
			else

			#Get Columns of mpileup and determine population number
				declare -i NCOL=$(awk '{ print NF; exit }' $OUTDIR/inds/${file}.mpileup)
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
				awk "$FOUR" $OUTDIR/inds/${file}.mpileup > $OUTDIR/inds/${file}_RindSTATs.rin 
			fi
		}

		N=$THREADZ
		open_sem $N
			for file in $(cat $OUTDIR/pops/${OUTPOP}_files_for_pops.txt) ; do
				run_with_lock task $i
		done ; wait
		echo "ALERT: Individual mpileups created at $(date) "

	#Create individual sync files
		echo "ALERT: Individual syncs being created at $(date) "

		for file in $(cat $OUTDIR/pops/${OUTPOP}_files_for_pops.txt) ; do
				if [[ -f ${OUTDIR}/inds/${file}.sync ]] ; then
					echo "ALERT: ${file}.sync exists, skipping this file "
				else
					syncin=$OUTDIR/inds/${file}.mpileup
					nice -n 19 java -ea -${KMEM} -Djava.io.tmpdir=`pwd`/tmp -jar ${POOL2}/mpileup2sync.jar --min-qual 1 --fastq-type sanger --input ${syncin} --output ${OUTDIR}/inds/${file}.sync --threads $THREADZ 
				fi
		done
		echo "ALERT: mpileup2sync.jar finished individual syncs at $(date) "

	#REQUIRES R: Get stats on individual files; these can be used to filter positions
		echo "ALERT: Rscript called to calculate population stats at $(date) "
		for file in $(cat $OUTDIR/pops/${OUTPOP}_files_for_pops.txt) ; do
				if [[ -f $OUTDIR/inds/${file}_snp_stats.txt ]] ; then
					echo "ALERT:${file}_snp_stats.txt  exists, skipping this analysis "
				else				
					rin=$OUTDIR/inds/${file}_RindSTATs.rin
					rout=$OUTDIR/inds/
					Rscript $BASEDIR/rscripts/r_ind_stats.R $rin $rout
				fi
		done
		echo "ALERT: Rscript stats done at $(date) "
		rm $OUTDIR/inds/*RindSTATs.rin

	#REQUIRES R: Standardize mpileups
		echo "ALERT: Rscript called to standardize syncs at $(date) "
		for file in $(cat $OUTDIR/pops/${OUTPOP}_files_for_pops.txt) ; do
				if [[ -f $OUTDIR/inds/${file}_norm.sync ]] ; then
					echo "ALERT:${file}_norm.sync exists, skipping this analysis "
				else				
					rin=$OUTDIR/inds/${file}.sync
					rout=$OUTDIR/inds/
					Rscript $BASEDIR/rscripts/r_standardize.R $rin $rout
				fi
		done
		echo "ALERT: Rscript standardization done at $(date) "

	#Combine black-listed polymorphic sites into one list of unique positions
				cat $OUTDIR/inds/*_poly_sites.txt | awk '{$3=$1+$2}1'| awk '!a[$3]++' | awk '{print $1,$2}' | sort -k1 -k2,2n > $OUTDIR/filters/${OUTPOP}_ind_poly_blacklist.txt
				if [[ -f $OUTDIR/inds/*_poly_sites.txt ]] ; then
					mv $OUTDIR/inds/*_poly_sites.txt $OUTDIR/filters/
				fi
		echo "ALERT: Normalized sync file order is: "
			ITER=1	
			for file in $OUTDIR/inds/*_norm.sync ; do
				namez=${file##*/}
  				echo ${ITER} = ${namez%*_norm.sync}
				ITER="$(($ITER + 1))"
			done	

	
	#Combine standardized mpileup into new sync
			if [[ -f $OUTDIR/${OUTPOP}_norm.sync ]] ; then
				echo "ALERT: Normalized sync file already exists; skipping"
			else
				for file in $(cat $OUTDIR/pops/${OUTPOP}_files_for_pops.txt| awk 'NR==1{print $1}') ; do
					paste <(awk '{print $1"\t"$2"\t"$3}' $OUTDIR/inds/${file}.mpileup) $OUTDIR/inds/*norm.sync >  $OUTDIR/${OUTPOP}_norm_temp.sync
				done
			fi

	#Filter indels from new sync by gtffile 
			if [[ -f ${OUTDIR}/${OUTPOP}_norm.sync ]] ; then
				echo "ALERT: Normalized indel file already exists; skipping"
			else
				perl ${POOL2}/indel_filtering/filter-sync-by-gtf.pl --input  $OUTDIR/${OUTPOP}_norm_temp.sync  \
						--gtf ${OUTDIR}/${OUTPOP}.gtffile --output ${OUTDIR}/${OUTPOP}_norm.sync
					declare -i before=$(wc -l < $OUTDIR/${OUTPOP}_norm_temp.sync)
					declare -i after=$(wc -l < ${OUTDIR}/${OUTPOP}_norm.sync)
					declare -i lost="$(($before - $after))"
				lostp=$((100-100*$after/$before))
					echo "ALERT: With an indel window of $INWIN bp you lost $lost SNPs or $lostp %" 
					rm $OUTDIR/${OUTPOP}_norm_temp.sync
			fi
fi

echo  "ALERT: Alignment and data creation step finished at $(date) "
echo "ALERT: Creating allele frequency tables in R"

#REQUIRES R: Frequency and MAF
		echo "ALERT: Rscript called to calculate allele frequencies at $(date) "
				if [[ -f ${OUTDIR}/${OUTPOP}_full.fz ]] ; then
					echo "ALERT: Frequency file ${OUTPOP}_full.fz  exists, skipping this analysis "
				else				
					rin=${OUTDIR}/${OUTPOP}.sync
					rout=$OUTDIR/
					rout2=$OUTDIR/filters/
					Rscript $BASEDIR/rscripts/r_frequency.R $rin $rout $rout2 $MAF
					echo "ALERT: Frequency file ${OUTPOP}_full.fz and its counterparts created at $(date) "
				fi
				
				if [[ ! -s ${OUTDIR}/${OUTPOP}.fz  ]]  ; then
					echo "ERROR: Frequency table is empty. Check for R errors"
					exit
				fi

				if [[ "$INDCONT" =~(on)$ ]] ; then

					if [[ -f ${OUTDIR}/${OUTPOP}_norm_full.fz ]] ; then
						echo "ALERT: Standardized Frequency file ${OUTPOP}_norm_full.fz exists, skipping this analysis "
					else				
						rin=${OUTDIR}/${OUTPOP}_norm.sync
						rout=$OUTDIR/
						rout2=$OUTDIR/filters/
						Rscript $BASEDIR/rscripts/r_frequency.R $rin $rout $rout2 $MAF $INDCONT
						if [[ ! -s ${OUTDIR}/${OUTPOP}_norm.fz  ]]  ; then
							echo "ERROR: Normalized frequency table is empty. You likely ran out of memory."
							echo "Consider reducing number of individuals or SNPs using higher filter thresholds."
							exit
						fi
						echo "ALERT: Frequency file ${OUTPOP}_norm_full.fz and its counterparts created at $(date) "
					
					fi
				fi

echo "ALERT: PPalign completed at $(date) "
			if [[ -f $OUTDIR/${OUTPOP}.gtffile.params ]] ; then
					rm $OUTDIR/${OUTPOP}.gtffile.params
			fi
			if [[ -f $OUTDIR/${OUTPOP}.sync.params ]] ; then
					rm $OUTDIR/${OUTPOP}.sync.params
			fi
			if [[ -f $OUTDIR/${OUTPOP}_norm.sync.params ]] ; then
					rm $OUTDIR/${OUTPOP}_norm.sync.params
			fi