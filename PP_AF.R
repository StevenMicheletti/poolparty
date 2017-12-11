#POOLSEQ SYNC to MAJOR ALLELE FREQUENCIES. V1. Steven Micheletti
#Takes a sync file and produces allele frequencies by population
#Designed for variant sites only, not entire genomes! (Unless they are < 10 mb)
#It is assumed at this point that you have applied a coverage filter

ppaf <- function(infile,
                      remove.uninformative = FALSE,
                      coverage.table = FALSE,
                      remove.poly = FALSE,
					  MAF=0.05,
					  maf.pos = FALSE) {

	#Dependencies
	require(plyr)
	require(dplyr)
	require(tidyr)
	require(fBasics)
	require(data.table)
	require(stringr)
	PF= 3

	#start timer
	stm <- proc.time()

	
#Some minor error checking
  if (MAF < 0 | MAF > 1) stop ("Minor allele frequency exceeds limits")

	print("Reading input file")
	system.time(gpops <- fread(infile,stringsAsFactors=FALSE))

	counts= gpops[, 4:(ncol(gpops))]
	gpops = gpops[, 1:3]  
	counts <- as.matrix(counts)

	activedir = getwd()
	tfile= paste0(activedir,"/","tempRf.tem")

	system.time(write.table(counts, file= tfile, sep=':',
						row.names=FALSE, col.names=FALSE, quote=FALSE))

	system.time(counts <- fread(tfile,stringsAsFactors=FALSE, sep=':'))

		file.remove(tfile)

	#Create rownames for freq file
	gpops3 <- unite(gpops, V1, c(V1,V2), sep="_", remove=FALSE)

	npops = ncol(counts)/6

		if (remove.poly == TRUE) {
		    print ("Removing sites with more than 3 alleles")
			counts$A <-rowSums(counts !=0)
			AL <- npops * PF
			rmrow <- nrow(counts[ which(counts$A> AL),])
			gpops3$A <- counts$A
			counts <- counts[ which(counts$A < AL + 1),]
			gpops3 <- gpops3[ which(gpops3$A < AL + 1),]
			counts$A=NULL
			gpops3$A=NULL
			print(c(rmrow,"loci removed for having more than 2 alelles"))
			} 

	#summary info from file
	totcol = ncol(counts)
	thetseq = seq(6,totcol,6)
	thebseq = seq(1,totcol,6)
	L = nrow(counts)


	#Break up Genotypes for each pop into separate lists
		Var <- list()
			for (i in (1:npops)) {
				Var[[i]] <-counts[, thebseq[i]:thetseq[i]]
		}

		Dcov <- data.frame(matrix(, nrow = L, ncol = npops))
			for (i in (1:length(Var))) {
			Dcov[[i]] <- rowSums(Var[[i]])
		}
		colnames(Dcov) <- c(1:npops)


#Apply Base Column Names
		Var<-llply(Var, function(x) {
		colnames(x) <- paste(c('A','T','C','G', 'N', 'del') )
		return(x) })

		MaxC <- list()
		Sumz<- llply(Var, function(x) rowSums(x), .progress='text')

		for (i in (1: npops)) {
			MaxC[[i]] <- colnames(Var[[i]])[max.col(Var[[i]],ties.method="first")]
			}

	print("Calculating frequencies... this will take a bit")
	Maxz<- llply(Var, function(x) rowMaxs(x), .progress='text')
	Freqz <- mapply("/",Maxz,Sumz,SIMPLIFY = FALSE)  


			for (i in (1: npops)) {
				MaxC[[i]] <- colnames(Var[[i]])[max.col(Var[[i]],ties.method="first")]
			}

	MaxC2 <- as.data.frame(MaxC)
	colnames(MaxC2) <- c(1:npops)
	MaxC3 <- as.matrix(MaxC2)

	Freqz2 <- as.data.frame(Freqz)
	colnames(Freqz2) <- c(1:npops)
	Freqz3 <- as.matrix(Freqz2)
	rownames(Freqz3) <- gpops3$V1

	PoStand <- as.data.frame(MaxC[[1]])
	colnames(PoStand) <- 'V1'
	PoStand3 <- as.matrix(PoStand)

	idx <- (!apply(MaxC3, 2, function(x) x == PoStand3))
	Freqz3[idx] <- 1- Freqz3[idx]
	Freqz3 <- na.omit(Freqz3)

	print("Frequency calculation done")

	#Write frequency table out
		if (remove.uninformative == TRUE) {
			print( "Removing positions with identical frequencies")
			d <- i+1
			DUpz <- as.matrix(rowStdevs(Freqz3))
			Freqz3 <- cbind(Freqz3, DUpz)
			Freqz3 <- subset(Freqz3, Freqz3[,d] != 0)
			Freqz3 <- (Freqz3[, -d])
			print( "Uniformative positions removed")
		}
		
		
	#Perform MAF filter
		mAx <- rowMaxs(Freqz3)
		mIn<- rowMins(Freqz3)
		rAnge <- as.vector(mAx - mIn)
		rAnge <- as.matrix(rAnge)
		rownames(rAnge) <- c()
		
		Freqz3 <- cbind (Freqz3, rAnge)

		
		if (maf.pos == TRUE) {
			print( "Getting positions of loci that do not meet MAF threshold")
			nfz <- ncol(Freqz3)
			alposs <- Freqz3[ which(Freqz3[,nfz] > MAF),]
			alposs2 <- as.data.frame(rownames(alposs))
			write.table(alposs2, file =tfile, sep='\t', quote=FALSE, col.names= FALSE, row.names = FALSE)
			gpout = c(infile,".affails")
			gpout = paste(gpout, collapse="")
			file.rename(tfile, gpout)
			print( "Position table created")
			}

			
			
nfz <- ncol(Freqz3)
FreqzX <- Freqz3[ which(Freqz3[,nfz] < MAF),]
Freqz3 <- Freqz3[ which(Freqz3[,nfz] > MAF),]
alert1 <- paste0((nrow(FreqzX)), " SNPs failed minor allele frequency filter")
print (alert1)

nfz <- ncol(Freqz3)
Freqz3 <- Freqz3[,-nfz]



	
	write.table(Freqz3, file= tfile, sep='\t',row.names= TRUE, quote=FALSE, col.names= TRUE)

	gpout = c(infile,".freq")
	gpout = paste(gpout, collapse="")
	file.rename(tfile, gpout)


		if (coverage.table == TRUE) {
			print("Creating coverage table")
			rownames(Dcov) <- gpops3$V1
			keepz <- as.vector(rownames(Freqz3))
			Dcov <- subset(Dcov, rownames(Dcov) %in% keepz)
			write.table(Dcov, file= tfile, sep='\t',row.names= TRUE, quote=FALSE, col.names= TRUE)
			gpout = c(infile,".cover")
			gpout = paste(gpout, collapse="")
			file.rename(tfile, gpout)
			print( "Coverage table created")
		}

	print ("Done! Look for '.freq' table")
proc.time() - stm
}
