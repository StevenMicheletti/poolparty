#!bin/R

suppressMessages(require(data.table))
suppressMessages(require(matrixStats))

#Get info for bash script
args <- commandArgs()
name <- args[6]
outdir <- args[7]
MAF <- as.numeric(args[8])

alertname <- sub('.*/', '', name )
alertname2 <- gsub( ".sync", "", alertname)

#Read file in
infile = fread(file=name, header=T, showProgress=FALSE)
infile = infile[complete.cases(infile), ]
freqz = as.matrix(infile[,4:ncol(infile)])

#Get MAF
	mAx <- rowMaxs(freqz, na.rm = TRUE)
	mIn <- rowMins(freqz, na.rm = TRUE)
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

#Subset
MAFFAIL <- as.data.frame(cbind (infile[,1:2],freqz, rAnge, mAx, mIn))
MAFFAIL2 <- subset(MAFFAIL[,1:2], ((MAFFAIL$MAX >= (1-MAF) | MAFFAIL$MIN <= MAF) & (MAFFAIL$MAFT <= MAF)) )

#Determine loss
alert1 <- paste0("R ALERT: Comparison-specific MAF filters (<", MAF,") removed ", nrow(MAFFAIL2), " SNPs")
print(alert1)

#Write out MAF blacklist
outname1 <- paste0(outdir, "/", alertname2, "_mafBL")
write.table(MAFFAIL2, file=outname1,
                row.names=FALSE, col.names=FALSE, quote=FALSE)