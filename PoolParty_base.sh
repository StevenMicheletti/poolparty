#!/bin/bash

#PoolParty_base_V1.0
#By Steven Micheletti
#1/3/2018 - remove echo from trim
#12/15/17

#Must have configuration file in same directory
source $PWD/PoolParty_base.config

declare -i THREADZ
declare -i QUAL
declare -i MINLENGTH
declare -i INWIN
declare -i MAPQ

#Copy configuration file so you know what you did later
cp $PWD/PoolParty_base.config ${RUNDIR}/${OUTPOP}.config


# MULTICOMP (or not)
if [[ "$MULTICORE" =~(on)$ ]] ; then
	echo "Multi-core is on, be weary of CPU/RAM usage"
else
 echo "Multi-core is off, will not parallel certain tasks"
fi

####-----------------------------------------RUNNING-----------------------------------------------#####
###At 30 threads, expect this whole process to take about 1 day per full Illumina Lane (at 400 million reads each)####

			startt=`date +%s`
			echo "Starting Quality Sequence Check"

#Check for Samplelist file, quit if it does not exist
	if [[ ! -f $INDIR/"samplelist.txt" ]] ; then
 	 echo 'File "samplelist.txt" is not there, aborting.'
   	 exit
		fi
		
#Check for index file, quit if it does not exist
	if [[ ! -f ${GENOME}.fai ]] ; then
 	 echo "${GENOME}.fai index file is missing, aborting"
   	 exit
		fi


#Sort sample list to get mate pairs in correct order
	sort -f $INDIR/"samplelist.txt" > $INDIR/"samplelist.tmp" && mv $INDIR/"samplelist.tmp" $INDIR/"samplelist.txt"
	FILES=$(cat $INDIR/"samplelist.txt")

#Make the output directory if it does not exist 
	if [ ! -d "$OUTDIR/" ]; then
		mkdir $OUTDIR/
	fi

#Check mate pairs; fail if order is incorrect
	cut -d_ -f1  $INDIR/"samplelist.txt" >  $OUTDIR/${OUTPOP}_prefixes.txt

	CHECK=$(cat $OUTDIR/${OUTPOP}_prefixes.txt)
	CHECK2=($CHECK)

	for (( i=0; i<${#CHECK2[@]} ; i+=2 )) ; do
		echo  completed pairs "${CHECK2[i+1]}" and "${CHECK2[i]}"
		if [ "${CHECK2[i+1]}" == "${CHECK2[i]}" ]; then
		echo "Correct Mates, proceeding"
	else
		echo "Incorrect samples, check samplelist ya dummy!"
    		exit
	fi
	done

#Create prefix file; this indicates the populations used and the order they are incorporated into population files
	cut -d_ -f1 $INDIR/"samplelist.txt" | awk '!seen[$0]++' | sed '/^$/d'  >  $OUTDIR/${OUTPOP}_prefixes.txt
	printf "using libraries:\n$(cat $OUTDIR/${OUTPOP}_prefixes.txt)\n"
	
#Get anchored chromosome lengths in bp. Prints a copy to the rundir
	echo "Creating genome anchored index..."
	HNAM=$(cut -f1 ${GENOME}.fai |  sed 's/[^a-Z]*//g'  |  awk '!x[$0]++')
	HNAM2=$(echo $HNAM)
		echo "$HNAM2 are your chromosome and/or scaffold headings"
	cut -f1,2 ${GENOME}.fai  | grep -v "^${SCAHEAD}"  >  $OUTDIR/${OUTPOP}_CHRbp1.txt
	awk '{print $1}' $OUTDIR/${OUTPOP}_CHRbp1.txt  | awk  '$2="1"' | awk '{gsub(" ","\t",$0); print;}' > $OUTDIR/${OUTPOP}_CHRbp2.txt
	cat $OUTDIR/${OUTPOP}_CHRbp1.txt $OUTDIR/${OUTPOP}_CHRbp2.txt > $OUTDIR/${OUTPOP}_CHRbp.txt
	cp  $OUTDIR/${OUTPOP}_CHRbp.txt $RUNDIR/CHRbp.txt &
	rm $OUTDIR/${OUTPOP}_CHRbp1.txt ; rm $OUTDIR/${OUTPOP}_CHRbp2.txt
	
#Trim by quality score 
	if [ ! -d "$OUTDIR/trimmed" ]; then
		mkdir $OUTDIR/trimmed
	fi
	
	array=($FILES)
	echo trimming  ${#array[@]} files 
	for (( i=0; i<${#array[@]} ; i+=2 )) ; do
		b=${array[i]%%_*} 
		if [[ -f $OUTDIR/trimmed/$b.trim ]] ; then
			echo 'trim file "$OUTDIR/trimmed/$b.trim" exists, aborting'
			exit
		fi
		
		if [[ "$MULTICORE" =~(on)$ ]] ; then
			nice -n 19 perl ${POPTRIM} --input1 $INDIR/${array[i]} --input2 $INDIR/${array[i+1]} --output $OUTDIR/trimmed/$b.trim --quality-threshold $QUAL --min-length $MINLENGTH --fastq-type $SCORETYPE &
			echo "Trimming pairs "${array[i+1]}" 1 and "${array[i]}" 2 " &
		else
			nice -n 19 perl ${POPTRIM} --input1 $INDIR/${array[i]} --input2 $INDIR/${array[i+1]} --output $OUTDIR/trimmed/$b.trim --quality-threshold $QUAL --min-length $MINLENGTH --fastq-type $SCORETYPE 
			echo "Trimming pairs "${array[i+1]}" 1 and "${array[i]}" 2 " 
		fi
		done;wait
			
# Quality report of trimmed files (FASTQC). This runs in the background as alignments are produced
	if [ ! -d "$OUTDIR/quality" ]; then
		mkdir $OUTDIR/quality
	fi
	
	if [[ "$MULTICORE" =~(on)$ ]] ; then
		for i in ${OUTDIR}/trimmed/*.{trim_1,trim_2}; do
			${FASTQC} -o ${OUTDIR}/quality $i &
		done &
	else
				for i in ${OUTDIR}/trimmed/*.{trim_1,trim_2}; do
			${FASTQC} -o ${OUTDIR}/quality $i 
		done &
	fi
	
#Alignment and duplicate removal (BWA, SAMBLASTER, SAMTOOLS)
#Discordant SAM files are printed as well for structural variant analyses
	#Make folder for BAMS
	
	if [ ! -d "$OUTDIR/BAM" ]; then
	mkdir $OUTDIR/BAM
	fi
	
	for i in $(cat $OUTDIR/${OUTPOP}_prefixes.txt); do
		if [[ -f ${OUTDIR}/BAM/$i.bam ]] ; then
   			echo 'BAM file "${OUTDIR}/BAM/$i.bam" exists, aborting'
			exit
		fi
		
		#Makes unique lane and pop name for each library
		a="@RG\tID:Lane"
		b="\tPL:ILLUMINA\tLB:test\tSM:"
		c=$a$i$b$i.1
		nice -n 19 ${BWA} mem -M -t $THREADZ -R $c \
		"${GENOME}" \
		"${OUTDIR}/trimmed/$i.trim_1" "${OUTDIR}/trimmed/$i.trim_2" | ${SAMBLASTER} -M -r -d "${OUTDIR}/BAM/$i.disc.sam" -s "${OUTDIR}/BAM/$i.split.sam"  | ${SAMTOOLS} view -Sb -q ${MAPQ} - > "${OUTDIR}/BAM/$i.bam"

	done; wait
		echo "FastQ Processing Sequence Completed"
		echo "Sorting and filtering out unused reads"

		#reports will have read alignment results for each bam file
		if [ ! -d "$OUTDIR/reports" ]; then
		mkdir $OUTDIR/reports
		fi

# Sorting BAM files and filtereing unpaired reads (PICARDTOOLS, SAMTOOLS)
	for i in $(cat $OUTDIR/${OUTPOP}_prefixes.txt); do
		nice -n 19 java -XX:ParallelGCThreads=$THREADZ -Xmx4g -Djava.io.tmpdir=`pwd`/tmp -jar ${PICARDTOOLS} SortSam I= ${OUTDIR}/BAM/$i.bam O= ${OUTDIR}/BAM/${i}_2.bam VALIDATION_STRINGENCY=SILENT SO=coordinate TMP_DIR=`pwd`/tmp 
	done; wait

	for i in $(cat $OUTDIR/${OUTPOP}_prefixes.txt); do
		if [[ "$MULTICORE" =~(on)$ ]] ; then
			nice -n 19 ${SAMTOOLS} flagstat ${OUTDIR}/BAM/${i}_2.bam > ${OUTDIR}/reports/${i}_${OUTPOP}_aln_report.txt &
		else
			nice -n 19 ${SAMTOOLS} flagstat ${OUTDIR}/BAM/${i}_2.bam > ${OUTDIR}/reports/${i}_${OUTPOP}_aln_report.txt 
		fi
		nice -n 19 ${SAMTOOLS} view -b -F 0x04 -f 0x02 ${OUTDIR}/BAM/${i}_2.bam  > ${OUTDIR}/BAM/${i}_3.bam 
	done; wait

		echo  "Data sorted and filtered"
		echo "Combining populations and performing poolseq analysis preparations"
		startp=`date +%s`
		
		
#Mpileup file: Combines BAM files in specified order (SAMTOOLS)

	COMBINE=$(sed 's,^,'$OUTDIR/BAM/', ; s/$/_3.bam/' $OUTDIR/${OUTPOP}_prefixes.txt)
	nice -n 19 ${SAMTOOLS} mpileup  -B $COMBINE | \
		awk '$4 >=1 && 7 >=1 && 10 >=1 && 13 >=1 && 16 >=1 && 19 >=1 && 21 >=1 && 24 >=1 && 27 >=1 && 30 >=1 && 33 >=1' | grep -Ev $'^\t|\t\t|\t$' > ${OUTDIR}/${OUTPOP}.mpileup
	

#Remove initial BAM files (comment out if you want to keep these)

	 for i in $(cat $OUTDIR/${OUTPOP}_prefixes.txt); do

		if [[ -s ${OUTDIR}/BAM/${i}.bam && -s ${OUTDIR}/BAM/${i}_2.bam ]] ; then
			echo Removing initial bam files
			rm ${OUTDIR}/BAM/${i}.bam
	else
		echo something was wrong with ${i}_.bam. Make sure your BAM files exist
	exit
	
	fi
done; wait

echo "Preparing Popoolation Files"

#In-del filtering 1: Identify indel regions arouund specified indel window (POPOOLATION2 identify-indel-regions.pl)
	echo "Hold up... making GTF file for indels"
	nice -n 19 perl ${INDELREG} --input ${OUTDIR}/${OUTPOP}.mpileup --indel-window $INWIN --output ${OUTDIR}/${OUTPOP}.gtffile --min-count 1

#Create Sync format for popoolation (POPOOLATION2 mpileup2sync.jar)
	echo "Hold up... making sync file"
	nice -n 19 java -ea -Xmx7g -jar ${MP2SYNC} --fastq-type sanger --min-qual $QUAL --input ${OUTDIR}/${OUTPOP}.mpileup --output ${OUTDIR}/${OUTPOP}.sync --threads $THREADZ 

#Indel filtering : Remove regions around indels in sync file (POPOOLATION2 filter-sync-by-gtf.pl)
	echo "Hold up... making removing indels from sync file"
	nice -n 19 perl ${FILTERSYNC} --input ${OUTDIR}/${OUTPOP}.sync --gtf ${OUTDIR}/${OUTPOP}.gtffile --output ${OUTDIR}/${OUTPOP}_indel.sync 

	if [[ -s ${OUTDIR}/${OUTPOP}.mpileup && -s ${OUTDIR}/${OUTPOP}.sync ]]
		then	
			rm ${OUTDIR}/${OUTPOP}.sync
		else
			echo Something is wrong with  ${OUTDIR}/${OUTPOP}.mpileup or ${OUTDIR}/${OUTPOP}.sync, exiting
			exit
		fi
			endp=`date +%s`
			echo  "Data prepared for analyses in $((endp-startp)) seconds, proceed to analyses... and... good luck..."

exit
