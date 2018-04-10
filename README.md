# :fish: PoolParty  :umbrella:

### A Pool-Seq Bioinformatic Pipeline (ver 0.76)

A BASH pipeline to align and analyze paired-end NGS data.

# Getting Started

Ensure that proper permissions are set to execute each package in the pipeline . 

 PoolParty is designed to be run on analysis servers. As such, memory and storage may be limiting factor for some systems depending on genome sizes, number of pools, number of SNPs, etc. In general, 128 gb of ram with at least 10 threads should suffice. 

 It is highly recommended to run the example files provided in the example directory before diving into large datasets.  

 PP_example.pdf under 'examples' contains additional detailed information on the pipeline.  


## Prerequisites

PoolParty is designed to be run on Linux (GNU) operating systems. Because it coordinates the execution of multiple packages there are number of dependencies that must be installed prior to running. With the use of diverse packages, the latest versions of Java, Perl, and R must be installed. The required packages for PoolParty are:

## Required package with version at inception 
- Burrows-Wheeler Aligner (BWA; 07.12) - http://bio-bwa.sourceforge.net/  
- Fastqc (0.11.7 ) - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/  
- samblaster (0.1.24) - https://github.com/GregoryFaust/samblaster  
- samtools (1.5) - http://www.htslib.org/download/  
- bcftools (1.5) - http://www.htslib.org/download/  
- Picard Tools (2.17.11) - http://broadinstitute.github.io/picard/  
- Popoolation2 (1.201) - https://sourceforge.net/p/popoolation2/wiki/Main/  
- BBMap (37.93) - https://sourceforge.net/projects/bbmap/  


### R packages used
If not already installed, PoolParty will attempt to automatically install required R packages. It is recommended to manually install packages beforehand:  

-PPalign: matrixStats, tidyr, stringr, data.table  
-PPstats: reshape, fBasics, ggplot2, RColorBrewer  
-PPanalyze: matrixStats, plyr, stringr, data.table, fBasics, ape, metap

## Input Files

PoolParty requires paired-end .fastq files (compressed or not) and a reference assembly genome (draft or complete). The reference genome must be indexed and prepared to properly perform alignments and sorting from various packages:

### Reference genome preparation:

Prepare the reference genome for bwa mem:
> $ bwa index -a bwtsw ref_genome.fasta

Index the reference genome for samtools/bcftools:
> $ samtools faidx ref_genome.fasta

Create dictionary for Picard Tools:
> $ -jar picard.jar CreateSequenceDictionary REFERENCE=ref_genome.fasta OUTPUT=ref_genome.fasta.dict



## Installing and Running The Pipeline 

Clone poolparty:

> $ git clone https://github.com/StevenMicheletti/poolparty/

There is no compiling or decompressing required. Place the folder into the desired directory:

> $ mv poolparty /usr/local/bin/poolparty

Creating a symbolic link will allow you to run each module from the command line :

> $ ln -s /usr/local/bin/poolparty/PPalign.sh /usr/local/bin/PPalign  
> $ ln -s /usr/local/bin/poolparty/PPanalyze.sh /usr/local/bin/PPanalyze  
> $ ln -s /usr/local/bin/poolparty/PPstats.sh /usr/local/bin/PPstats  

To run a PoolParty module in unix followed by a corresponding configuration file :

> $ PPalign align.config

For each run, it is highly recommended to run the pipeline into a log file :

> $ PPalign align.config &> run.log &


# PPalign

PPalign takes paired-end fastq files through the process of quality filtering and alignment to a reference genome.
The specified populations are combined into a VCF file, a mpileup file, a sync file, and independent bam files. Additional reports are generated as well.

## What is the alignment phase doing to your fastqs? 

If pooled data do not have individual barcodes:
- Quality trimming of fastq files and contaminant removal (BBduk from BBMap)  
- Producing quality summaries of trimmed fastqs (fastqc)  
- Aligning trimmed fastq files to genome assembly (bwa mem)  
- Removing PCR duplicates and producing discordant and split-end bam files (samblaster)  
- Filtering the aligned bam files and producing alignment stats (samtools)  
- Sorting the bam files (Picard)  
- Calling SNPs across all libraries (bcftools)
- Identifying in/del regions  (Popoolation2)  
- Creating sync file with in/del regions masked (Popoolation2)  
- MAF filtering and potential paralog identification (r_frequency.R)

If individual barcoded analyses the alignment phase also contains:

- Individual contribution stats (r_ind_stats.R)
- A standardized sync file that normalizes individuals' allelic contribution to each genomic position (r_standardize.R)

## Input files
In addition to fastq files, PPalign also requires an input file named 'samplelist.txt' which must be placed in the directory containing the fastqs. samplelist.txt is a list containing the file names (one per line) for all of the fastqs you want to include in the run. An example with 2 paired-end libraries:

>Pool1_R1.fq.gz	1  
>Pool1_R2.fq.gz	1  
>Pool2_R1.fq.gz	2  
>Pool2_R2.fq.gz	2  

Note that the 'library' or 'population' number must be specified after each file. This number must be an integer from 1-N. See the example samplelist.txt for more information.

The naming convention of the fastq files is essential. The unique ID identifying the library must occur before the first underscore and must match its paired-end mate. The number after the file designates the population or library that the file belongs to - this is particularly useful if individuals are barcoded or populations were sequenced on different lanes. 


## Editing the .config file

The configuration file contains working directory locations, run parameters and dependency locations. There are example config files provided. Create a new .config file for each run and place it in your working directory.

#### Directory and input explanation

- INDIR = (dir; required) the input directory that contains unzipped fastq files. 'samplelist.txt' must also be here  
- OUTDIR = (dir; required) the base directory where output files will be written to  
- OUTPOP = (string; required) the unique prefix name for your population output files  
- GENOME = (file; required) the location and name of the fasta genome file
- SCAHEAD = (string; optional) the prefix that identifies unanchored scaffolds in the genome assembly  

#### Run Parameters

- THREADZ=(integer; required) the number of threads to use when multi-threading and parallel processing is available  
- BQUAL=(integer; required) minimum PHRED base quality for trimming raw reads 
- MAPQ=(integer; required) minimum MAPQ score to retain an alignment 
- SNPQ=(integer; required) minimum bcftools snp QUAL score to retain a SNP
- MINLENGTH=(integer; required) minimum length a fastq read can be trimmed to before discarding
- INWIN=(integer; required) the indel window size (bp) for masking SNPS around an indel 
- MAF=(float; required) minimum global minor allele frequency to retain a SNP
- KMEM=(string; required) maximum memory allocation for java packages. Value will vary based on system 
- MINDP=(integer; required) minimum global coverage needed to retain a SNP

#### Run-types
- SPLITDISC=(off/on; required) if on, produces split-end and discordant sam files. Not recommended unless the goal is to look at structural variants  
- INDCONT=(off/on; required) if on, will analyze fastqs as if they are independent individuals. Individual stats and normalization. Note this a high memory process  
- QUALREPORT=(off/on; required) if on, will produce quality reports from fastqc for each fastq file 

#### Dependency Locations
Identify the location and names of the executables / scripts.  If you've made programs executable across from any directory you don't need to include the directory  

- BCFTOOLS (executable) = bcftools
- FASTQC (executable) = fastqc
- BWA (executable) = bwa
- SAMBLASTER (executable) = samblaster
- SAMTOOLS (executable) = samtools
- PICARDTOOLS (file) = picard.jar
- BBMAPDIR (directory) = location of the bbmap directory
- POOL2 (directory) = location of the Popoolation2 directory

## Output files and directories
Many files will be produced during the alignment phase. Ensure you have enough storage before executing.

- ##### OUTDIR/OUTPOP_CHRbp.txt  
  - Contains the end position and start position in basepairs for each  chromosome
  
- ##### OUTDIR/OUTPOP.mpileup  
  - Combined filtered bam files, variants only. Order of columns is the same as in the order pools are listed in OUTDIR/_names.txt
  
- ##### OUTDIR/OUTPOP_stats.mpileup  
  - Coverage for every combined library across all genomic positions of the ref assembly. Used to generate mapping stats

- ##### OUTDIR/OUTPOP.VCF 
  - Filtered variant call format of all libraries consisting of variants only

- ##### OUTDIR/OUTPOP_variants.txt  
  - List of positions of variants from VCF file with QUAL and DP
  
- ##### OUTDIR/OUTPOP.sync  
  - Sync format with indel regions masked. Used for downstream analyses  
  
- ##### OUTDIR/_.fz
  - Allele frequency table of SNPs. Full is table with missing data, complete is table without missing data  
  
- ##### OUTDIR/trimmed/
  - Quality trimmed fastq files written here with corresponding trim reports

- ##### OUTDIR/quality/  
  - fastqc quality summaries for trimmed reads  

- ##### OUTDIR/BAM/.bam
  - Aligned (bwa-mem, samblaster duplicate removal), coordinate-sorted (picard), and filtered (samtools unpaired read removal) bam  files for each library. Filtered bams are also combined by population/category assignment from the samplelist

- ##### OUTDIR/BAM/split.sam  
  - If specified, split-end alignments produced by samblaster. Required if running SV-analysis such as LUMPY

- ##### OUTDIR/BAM/disc.sam  
  - discordant alignments produced by samblaster. Required if running SV-analysis such as LUMPY  

- ##### OUTDIR/reports/
  - Read alignment reports based on the alignment of trimmed reads for each aligned bam file
  
- ##### OUTDIR/pops/
  - Files which specify which library belongs to which population. Also indicates the order of libraries in each population

- ##### OUTDIR/filters/
  - Files that specify coverage for each library/population and SNP and INDEL genomic locations. Blacklisted MAF and > 2 allele SNPs will be produced here as well.
  
- ##### OUTDIR/inds/
  - If individual analyses turned on, individual-based stats, sync files, and mpileup files will be produced here. Additionally, normalized files will appear in OUTDIR/

# PPanalyze

## Input

PPanalyze uses a freq and sync file generated by PPalign to perform basic comparative analyses.

## How is PPanalyze analyzing you data? 

- FST and/or sliding window FST (Popoolation2)  
- Fisher's exact test for allele frequency differences (Popoolation2)
- Neighbor-joining trees (r_structure.R, ape)
- SNP density (r_structure.R)

## Editing the .config file

The configuration file contains working directory locations, run parameters and dependency locations. Create a new .config file for each run and place it in your working directory.

#### Directory and input explanation
- POPS=(alphanumeric string; required) - populations/libraries you wish to analyze/compare to one another. If more than two populations, comparative analyses (such as FST) will be averaged across all comparisons. In some cases, populations may share a same trait of interest and should not be averaged. A comma (,) between populations means compare those populations, a colon (:) means ignore that comparison.
- PREFIX=(string; required) - unique prefix for output files
- COVFILE=(file; required) - file produced by PPalign which contains depth of coverage information for all populations in the analysis
- SYNC=(file; required) - sync file produced by PPalign
- FZFILE=(file; required) - allele frequency file produced by PPalign
- BLACKLIST=(file; optional) - optional list of loci to black list (CHR, POS). Usually this is an output file from PPalign
- OUTDIR=(dir; required) - location to output files

#### Types of analyses
- FST=(on/off; required) - perform FST on each SNP
- SLIDINGFST=(on/off; required) - perform sliding window FST on defined SNP window
- FET=(on/off; required) - perform Fisher's exact test on each SNP
- NJTREE=(on/off; required) - perform NJ-tree and SNP density either on windowed or single SNPs

#### Parameters
- MINCOV=(integer; required) - The minimum coverage at a SNP (across every population included in the analysis) for it to be retained for any analysis
- MAXCOV=(integer; required) - The maximum coverage at a SNP (across every population included in the analysis) for it to be retained for any analysis
- MAF=(float; required) - comparison-specific minimum minor allele frequency to retain SNPs based on
- FSTTYPE=(traditional/karlsson; required IF FST=on) - Perform traditional FST analysis or Karlsson (et al. 2007) method (Uses allele depth) 
- WINDOW=(integer; required IF SLIDINGFST=on OR NJTREE=on) - window size for windowed approaches such as FST and NJTREE. If set to 1, NJTREE will not perform windowed analyses.
- STEP=(integer; required IF SLIDIGNFST=on) - step size for sliding window FST
- NIND=(integer; required IF FST=on) - number of individuals in pool for FST correction
- BSTRAP=(integer; required IF NJTREE=on) - number of bootstraps to perform on neighbor joining tree
- AFILT=(float; required if NJTREE=on ) - crude filter to try to eliminate SNPs under selection. Filters out SNPs with allele frequencies greater than this value. If set to 1, will not perform allele frequency filtering.
- METHOD=("mean","sdmax","sdmin","random","rangemax","rangemin","first"; required if NJTREE=on AND WINDOW >1) - Method to choose SNPs within a window for NJTREE analyses.

#### Dependency location
- POOL2 (directory) = location of the Popoolation2 directory
Note that FET analyses required an additional perl module which can be installed with:
> $ cpan Text::NSP::Measures::2D::Fisher::twotailed


## Output files 

- ##### OUTDIR/PREFIX_raw*
  - Popoolation2's raw output for the given analysis

- ##### OUTDIR/PREFIX_analysis*.txt
  - Reformatted output for the given analysis. CHR,POS, ID, and P/FST 
  
- ##### OUTDIR/PREFIX.sync
  - Subset sync file for specified populations
  
- ##### OUTDIR/PREFIX.fz
  - Subset frequency table for specified populations
  
- ##### OUTDIR/PREFIX_density.txt
  - SNP density for specified window with median position in the window
  
- ##### OUTDIR/PREFIX_consensus/single.pdf
  - If specified, NJ tree for a single tree and consensus tree using Nei's genetic distance and specified number of bootstraps
  
# PPstats

PPstats uses a mpileup file to perform depth of coverage statistics. This is particularly useful to assess sequencing performance.

## How does PPstats generate statistics?

PPStats simply takes a mpileup with each population's depth of coverage for each genomic position and determines stats such as:  
  -Mean depth of coverage for each population  
  -Mean depth of coverage after minimum and maximum coverage filters are applied for each population  
  -Proportion of the reference genome that each population covers with sufficient coverage  
  -Proportion of each chromosome covered (checking for biased alignment)  

## Editing the .config file

- FAI==(file; required) - index file for the genome. Contains chromosome/scaffold lengths. Produced by samtools faidx
- MPILEUP=(file; required) - stats mpileup generated by PPalign
- OUTDIR=(dir; required) - output directory
- OUTFILE=(string; required) - Prefix name of the stats output file (with no directory)
- SCAFP=(string; required) - prefix name for scaffolds (will be indicated in .fai file). This is essential for plotting chromosome vs scaffold results. If no chromosomes are anchored, chromosomal analyses will be bypassed.
- THREADZ=(integer; required) - number of threads to use. Memory usage is low for these analyses.
- MINCOV=(integer; required) - minimum desired coverage to retain a genomic position
- MAXCOV=(integer; required) - maximum desired coverage to retain a genomic position

## Output files 

- ##### OUTDIR/_summary.txt
  - Ref assembly stats
  
- ##### OUTDIR/_prop_cov.pdf
  - Proportion of genome covered after specified min and max depth of coverage
  
- ##### OUTDIR/_prop_at_covs.pdf
  - Proportion of genome covered at minimum depth increments
  
- ##### OUTDIR/_mean_filt_cov.pdf
  - Mean depth of coverage after filtering out reads by min and max depth of coverage
  
- ##### OUTDIR/_mean_coverage.pdf
  - Mean depth of coverage of all mapped reads
  
- ##### OUTDIR/_chr_prop_mean.pdf
  - Mean proportion of each chromosome covered by mapped reads after filtering out reads by min and max depth of coverage. Similar plots will produced for each library/population
 

# PPanalyze

Multiple utility scripts come with PoolParty for plotting and additional analyses. Each of these can be run from the terminal without the need for a configuration file. Typing -h after calling the script will bring up help options for each.

## PPmanhat
 - Plots results such as FST, SFST, FET, density in a 4-column format (CHR,BP,SNPID,FST)
## PPrunls
 - Runs local score (Fariello et al.) on p-values from SNPs
## PPruncmh
 - Runs a Cochran–Mantel–Haenszel (CMH) on .sync file
## PPsubset
 - Subsets FST,SFST, or FET file by specific libraries are additional coverage thresholds 

# Troubleshooting

PoolParty is new so users may encounter bugs. However, there are common issues that can be avoided  :

1) Permissions: Proper permissions are not only needed to run the PoolParty modules, but also all dependencies. Ensure that your user account has permissions to execute programs and write to the specified output directories.
2) Memory: With increased data comes increased memory usage. If java programs encounter a memory error they will usually spit out an interpretable error. Tune the java memory parameter accordingly.
3) Storage: Large temporary files can fill up smaller hard drives fast. Storage issues generally will have to be resolved with hardware. 
4) Compatibility: PoolParty is POSIX compliant but incompatibilities with specific Linux distributions can still be encountered. Specific formatting of output drives can cause issues, especially if piping is not supported on these drives (mkfifo). Errors associating with drives may require reformatting or diverting output files to a different drive. 

If an issue does not fall within this category, post the error message and explanation to the PoolParty GitHub page. Also, don't forget to check out the example file for more details. 

## Stopping PoolParty

Perhaps you included the wrong samples or need to add an additional information to a PoolParty run that is currently underway. Since many modules run background processes, you will have to kill the entire script in this fashion: 

Enter the module that is running and determine processes it is associated with:  

> $ ps -aux | grep PPalign 

Kill all script processes:  

> $ killall PPalign

You may need to kill any additional lingering processes

> $ kill PID
