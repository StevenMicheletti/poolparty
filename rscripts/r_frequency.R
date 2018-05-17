#!bin/R
#load packages
suppressMessages(require(tidyr))
suppressMessages(require(data.table))
suppressMessages(require(stringr))
suppressMessages(require(matrixStats))

#Get info for bash script
args <- commandArgs()
name <- args[6]
outdir <- args[7]
outdirF <- args[8]
MAF<- as.numeric(args[9])
#Names for output
alertname <- sub('.*/', '', name )
alertname2 <- gsub( ".sync", "", alertname)
#Read file 
alert1<-paste0("R ALERT: Reading file ", alertname, " into R")
  #Reading input .sync file
  gpops <- fread(file=name,stringsAsFactors=FALSE, showProgress=FALSE)
  #Break into genotypes and genomic positions heading
  counts= gpops[, 4:(ncol(gpops))]
  gpops = gpops[, 1:3] 
  colnames(gpops) <- c("Chr", "Pos", "Ref")
  #Convert to matrix for faster write
  counts <- as.matrix(counts)
  #Get writing dir, make temp file
  activedir = getwd()
  tfile= paste0(outdir,"/","tempRf.tem")
  #Write temp file of allele counts, this allows to break up colon-sep alleles
  write.table(counts, file= tfile, sep=':',
     row.names=FALSE, col.names=FALSE, quote=FALSE)
  #Read counts back in separate columns
  counts <- fread(tfile,stringsAsFactors=FALSE, sep=':', showProgress=FALSE)
  #Remove temporary file
  invisible(file.remove(tfile))
  print("R ALERT: Formatting file and calculating summary stats")
  #Get number of pops in this analysis (sync has 6 columns per pop)
  npops = ncol(counts)/6
  #Get summary info from file
  totcol = ncol(counts)
  thetseq = seq(6,totcol,6)
  thebseq = seq(1,totcol,6)
  L = nrow(counts)
  #Break up Genotypes for each individual into separate lists
  Var <- list()
    for (i in (1:npops)) {
      Var[[i]] <-counts[, thebseq[i]:thetseq[i]]
      gc()
    }
  #Get coverage at position ( skip if normalized file)
    if (is.na(args[10])) {
        print("R ALERT: Calculating depth of coverage at each position")
        Dcov <- data.frame(matrix(, nrow = L, ncol = npops))
        for (i in (1:length(Var))) {
           Dcov[[i]] <- rowSums(Var[[i]])
           gc()
        }
        #Determine N and additional indel sites
        indel.N = list()
        for (i in (1: npops)) {
             indel.N[[i]] <- as.matrix(rowSums(Var[[i]][,5:6]))
        }
        indel.N <- as.data.frame(indel.N)
        colnames(indel.N) <- c(1:npops)
        colnames(Dcov) <- c(1:npops)
        Dcov$Total <- rowSums(Dcov)
        Dcov <- cbind(gpops[,1:2], Dcov, indel.N)
        outname7 <- paste0(outdirF, alertname2, "_coverage.txt")
        write.table(Dcov, file=outname7,
               row.names=FALSE, col.names=TRUE, quote=FALSE)  
    }

#Frequency and paralog detection 
alertname <- sub('.*/', '', alertname )
  #Apply columns names as bases
  Var<-lapply(Var, function(x) {
     colnames(x) <- paste(c('A','T','C','G', 'N', 'del') )
  return(x) })
  #Make a numerical matrix for mode calculation
  VarN<-lapply(Var, function(x) {
    colnames(x) <- paste(c('1','2','3','4', '5', '6') )
    return(x) })
  print("R ALERT: Noting potential paralogs (>3 alleles per position)")
  #Determine polymorphic sites (>3 alleles). Generally memory heavy...can be commented out 
    MCNT <- list()
    MVar <- list()
    for (i in (1: npops)) {
     MVar[[i]] <- as.matrix(Var[[i]][,1:4])
      MCNT[[i]] <- rowSums(MVar[[i]] >0 )
      MCNT[[i]] <- as.data.frame(MCNT[[i]])
      MCNT[[i]] <-rapply(MCNT[[i]], function(x) ifelse(x > 2,1,0), how = "replace")
      gc()
    }
  PolY <-  Reduce("+", MCNT)
  rm (MCNT) ; invisible(gc())
  PolY <- cbind(gpops[,1:2], PolY)
    #Polymorphic in at least one populations
      PolY1<- PolY[,1:2][ PolY[[3]] > 0]
    #Polymorphic in at least half of the populations
      PolY2<- PolY[,1:2][PolY[[3]] >= npops*0.5]
    #Polymorphic in all populations
      PolY3<- PolY[,1:2][PolY[[3]] >= npops]
  outname4 <- paste0(outdirF, alertname2, "_poly_one.txt")
  outname5 <- paste0(outdirF, alertname2, "_poly_half.txt")
  outname6 <- paste0(outdirF, alertname2, "_poly_all.txt")
    write.table(PolY1, file=outname4,
                row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(PolY2, file=outname5,
                row.names=FALSE, col.names=FALSE, quote=FALSE) 
    write.table(PolY3, file=outname6,
                row.names=FALSE, col.names=FALSE, quote=FALSE)
    rm(PolY) ; invisible(gc())
    rm(PolY1) ; invisible(gc())
    rm(PolY2) ; invisible(gc())
    rm(PolY3) ; invisible(gc())

  #Get the major base at each position in each population
  #If a lot of missing data at a position, this will default to the ref base
  MaxC <- list()
    for (i in (1: npops)) {
       MaxC[[i]] <- colnames(Var[[i]])[max.col(Var[[i]],ties.method="first")]
       gc()
    }
  #Get the major base at each position in numerical set
  MaxN <- list()
    for (i in (1: npops)) {
       MaxN[[i]] <- colnames(VarN[[i]])[max.col(VarN[[i]],ties.method="first")]
       gc()
    }
  #Make nucleotide list into data frame, printing the primary allele for each population
    MaxC <- as.data.frame(MaxC)
    colnames(MaxC) <- c(1:npops)
    MaxC <- as.matrix(MaxC)
  #Make nucleotide list into data frame, printing the primary allele for each population
    MaxN <- as.data.frame(MaxN)
    colnames(MaxN) <- c(1:npops)
    MaxN <- as.matrix(MaxN)
    MaxN <- data.frame(apply(MaxN, 2, function(x) as.numeric(as.character(x))))
  #Get total coverage at positions 
  Sumz<- lapply(Var, function(x) rowSums(x))
  #Make into data matrix
  for (i in (1:npops)) {
    Var[[i]] <- as.matrix(Var[[i]])
  }
  #Calculate frequency of major allele
  print("R ALERT: Calculating Allele Frequencies")
  Maxz<- lapply(Var, function(x) rowMaxs(x))
  Freqz <- mapply("/",Maxz,Sumz,SIMPLIFY = FALSE)  

  #Merge into frequency table of major allele
  Freqz <- as.data.frame(Freqz)
  colnames(Freqz) <- c(1:npops)
  Freqz <- as.matrix(Freqz)
  
  #Reference allele based on assembly
  PoStand <- as.data.frame(gpops$Ref)
  colnames(PoStand) <- 'V1'
  PoStand <- as.matrix(PoStand)
  
  #Make reference alleles based on most seen allele in data set

  modefunc <- function(x){
    tabresult <- tabulate(x)
    themode <- which(tabresult == max(tabresult))
    if(sum(tabresult == max(tabresult))>1) themode <- themode[[1]]
    return(themode)
  }
  
  #Find most common allele, if >1, just print the first
  MaxN <- apply(MaxN, 1, modefunc)
  
  MaxN <- as.data.frame(chartr("123456", "ATCGND", MaxN))
  ComPos <- as.matrix(MaxN)
  colnames(ComPos) <- "A1"
  
  #Using reference allele,create comparisons index based on matching alleles
    #Make index, TRUE if the reference major allele is different in population x
      #Index based on most common allele
      idx <- (!apply(MaxC, 2, function(x) x == ComPos))
    
  #Determine alternative allele
     #If major allele isn't same between populations, make it the alternative allele
     Freqz[idx] <- 1- Freqz[idx]
     Freqz<- round(Freqz, 3)
     
    print("R ALERT: Determining Positions that fail MAF")
    #Get minor allele frequency failing positions
    mAx <- rowMaxs(Freqz, na.rm = TRUE)
    mIn<- rowMins(Freqz, na.rm = TRUE)
    rAnge <- as.vector(mAx - mIn)
    rAnge <- as.matrix(rAnge)
    mIn <- as.matrix(mIn)
    mAx <- as.matrix(mAx)
    rownames(rAnge) <- c()
    colnames(rAnge) <- "MAFT"
    rownames(mAx) <- c()
    colnames(mAx) <- "MAX"
    rownames(mIn) <- c()
    colnames(mIn) <- "MIN"
    
    MAFFAIL <- cbind (gpops[,1:2],Freqz, rAnge, mAx, mIn)
    MAFFAIL <- as.data.frame(MAFFAIL)
    MAFFAIL <- subset(MAFFAIL[,1:2], ((MAFFAIL$MAX >= (1-MAF) | MAFFAIL$MIN <= MAF) & (MAFFAIL$MAFT <= MAF)) )

    proP<- round((nrow(MAFFAIL) / nrow(Freqz) * 100),2)
    paste0("R ALERT: ", nrow(MAFFAIL), " SNPS (", proP, "%) do not pass additional population MAF threshold of ", MAF, " and have been noted")

#Create Table
    print("R ALERT: Writing output files")

    Freqz <- cbind(gpops[,1:3],ComPos,Freqz)

# Make a frequency table with non missing data
outname1 <- paste0(outdir, alertname2, ".fz")
outname3 <- paste0(outdirF, alertname2, "_MAF_fail.txt")
outname7 <- paste0(outdirF, alertname2, "_indelN.txt")

    #Write temp file of allele counts
    write.table(Freqz, file=outname1,
                row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(MAFFAIL, file=outname3,
                row.names=FALSE, col.names=FALSE, quote=FALSE)