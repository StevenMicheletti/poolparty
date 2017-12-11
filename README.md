# :fish: PoolParty :umbrella:

### A Pool-Seq Bioinformatic Pipeline

PoolParty takes raw paired-end fastq files, filters them, and formats them appropriately for further analysis while providing some stats along the way. Associated R scripts perform pair-wise analyses

## Prerequisites

PoolParty is designed to be run on Linux operating systems and primarily uses Unix tools. Because it coordinates the execution of multiple packages there are number of dependencies that must be installed prior to running. With the use of diverse packages, the latest versions of Java, Perl, and Python must be installed. The required packages for PoolParty are:

- Burrows-Wheeler Aligner (BWA) - http://bio-bwa.sourceforge.net/  
- FASTQC - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/  
- pyfaidx - https://pypi.python.org/pypi/pyfaidx  
- samblaster - https://github.com/GregoryFaust/samblaster  
- samtools - http://www.htslib.org/download/  
- Picard Tools - http://broadinstitute.github.io/picard/  
- Popoolation - https://sourceforge.net/projects/popoolation/  
- Popoolation2- https://sourceforge.net/p/popoolation2/wiki/Main/  

## Input Files

PoolParty requires unzipped .fastq files and a reference assembly genome. The reference genome must be thoroughly indexed to properly perform alignments and sorting:

> $ bwa index -a bwtsw  
> $ samtools faidx  
> $ picard.jar CreateSequenceDictionary  

PoolParty also requires an input file named 'samplelist.txt' which must be placed in the directory containing the fastqs. samplelist.txt is a list containing the file names (one per line) for all of the fastqs you want to include in the run. An example with 2 paired-end libraries:

>Pool1_R1.fastq  
>Pool1_R2.fastq  
>Pool2_R1.fastq  
>Pool2_R2.fastq  

The naming convention of the fastq files is essential. The unique ID identifying the library must occur before the first underscore and must match its paired-end mate. 

## Running The Pipeline 

PoolParty_base.sh and PoolParty_base.config must be in the same directory. Edit PoolParty_base.config and fill-in dependency and directory locations as well as parameter values. Simple execution of the alignment phase of the pipeline:

> $ ./PoolParty_base.sh

However, it is recommended to save the log as it will contain some valuable statistics...

> $ ./PoolParty_base.sh &> run.log &

or run in the background 

> $ nohup nice -n 19 ./PoolParty_base.sh &> run.log &

## Processes

What is the alignment phase doing to your fastqs? 

- Quality trimming fastq files (Popoolation)  
- Producing quality summaries of trimmed fastqs (fastqc)  
- Making a genome index of the reference genome (pyfaidx)  
- Aligning trimmed fastq files to genome assembly (bwa mem)  
- Removing PCR duplicates and producing discordant and split-end BAMs (samblaster)  
- Filtering the aligned BAM files and producing alignment stats (samtools)  
- Sorting the BAM files (Picard)  
- Combining filtered BAM files into mpileup (samtools)  
- Identifying in/del regions in mpileup (Popoolation)  
- Creating sync file with in/del regions masked (Popoolation)  

## Editing the .config file

The configuration file contains working directory locations, run paramteres and dependency locations.

#### Directories and inputs

- INDIR = (dir) the input directory that contains unzipped fastq files. 'samplelist.txt' must also be here
- OUTDIR = (dir) the base directory that contains output files  
- RUNDIR =(dir) the directory where you are running the bash script from  
- OUTPOP = (string) the unique prefix name for your population output files  
- GENOME = (file) the location and name of the fasta genome file
- SCAHEAD = (string) the prefix that identifies unanchored scaffolds in the genoem assembly. Usually "scaff" or the like  

#### Run Parameters

- THREADZ=(integer) the number of threads to use when multi-threading is possible  
- QUAL=(integer) minimum PHRED base quality 
- MINLENGTH=(integer) minimum length a fastq read can be trimmed to before throwing it out 
- INWIN=(integer) the indel window size (bp) for masking SNPS around an indel 
- MAPQ=(integer) minimum MAPQ score to retain an alignment 
- MULTICORE= ("on" or "off"). Turn on at your own risk. This will speed up your run by running certain scripts in parallel:
   trimming of fastq files, FASTQC quality scores, and samtools alignment stats. It roughly uses one core per fastq file. This may be      ideal for runs containing small numbers of pools, but may eat up RAM/CPU if you are running > 10 pools at once.  
   
#### Dependency Locations
Identify the location and names of the executables / scripts.  If you've made programs executable across the whole system you don't need to include the directory.

- FAIDX (file) faidx  
- FASTQC (file) fastqc  
- BWA (file) bwa  
- SAMBLASTER (file) samblaster  
- SAMTOOLS (file) samtools  
- PICARDTOOLS (file) picard.jar (Picard Tools, Java)  
- POPTRIM (file) trim-fastq.pl (Part of Popoolation, Perl)  
- INDELREG= (file) identify-indel-regions.pl (Part of Popoolation, Perl)  
- MP2SYNC (file) mpileup2sync.jar (Part of Popoolation, Java)  
- FILTERSYNC (file) filter-sync-by-gtf.pl (Part of Popoolation, Perl)  

## Output files and directories
Many files will be produced during the alignment phase. Ensure you have enough storage before executing.

- ##### OUTDIR/trimmed/trim_1  
  - Quality trimmed versions of the input fastq files. These are what get aligned to the genome assembly 

- ##### OUTDIR/OUTPOP_prefixes.txt  
  - The prefix names (Library ID) of your libraries for this run 

- ##### OUTDIR/OUTPOP_CHRbp.txt  
  - Contains the end position and start position (would should be 1) in basepairs for each anchored chromosome.

- ##### OUTDIR/quality/fastq  
  - fastqc quality summaries for trimmed reads  

- ##### OUTDIR/BAM/BAM_2  
  - Aligned, MAPQ filtered, and coordinate sorted BAM files  

- ##### OUTDIR/BAM/BAM_3  
  - Filtered BAM_2. All discordant, split-end, unpaired reads are filtered   

- ##### OUTDIR/BAM/split.sam  
  - split-end alignments produced by samblaster. Important if running SV-analysis such as LUMPY

- ##### OUTDIR/BAM/disc.sam  
  - discordant alignments produced by samblaster. Important if running SV-analysis such as LUMPY  

- ##### OUTDIR/OUTPOP.mpileup  
  - Combined BAM_3 files. Order of columns is the same as in the order pools are listed in the prefix file  

- ##### OUTDIR/OUTPOP_indel.sync  
  - Sync format with indel regions masked. Used for downstream analyses  

- ##### OUTDIR/reports/aln_report.txt
  - Read alignment reports based on the alignment of trimmmed reads
  
- ##### RUNDIR/.log
  - The run log not only contains run information, but trimming stats, duplicate stats, and alignment stats.
    
 ## Analysis scripts

- Coming soon

### PP_AF.R
- This R function (ppaf)  takes a .sync file and coverts it into an allele frequency table. It additionally will remove and produce a list of genomic positions failing minor allele frequency thresholds. If desired, a coverage table of the variant sites for each population/pool can also be produced.


