#!bin/R
suppressMessages(require(data.table))
suppressMessages(require(metap))

#Combine pval, transform back to -log10p

#remove sci notation
options(scipen=999)

args <- commandArgs()
name <- args[6]
outdir <- args[7]

	alertname <- sub('.*/', '', name )
	alertname2 <- sub('\\..*', '', alertname)

snps <- fread(name, stringsAsFactors=FALSE, header=F, showProgress=FALSE)
snps <- as.data.frame(snps)
snps [is.na(snps)] <- 0
snps[snps=="na"]<-0
snps <- apply(snps,2, function(x) as.numeric(as.character(x))) 


#Get p-val for each column
if (ncol(snps) > 1) {
	snps3 <- apply(snps,2, function(x) 10^-(x))
	snps4 <- as.matrix(apply(snps3,1, function(x) sumlog(x)[[3]])) 
}
if (ncol(snps) < 2) {
	snps3 <- apply(snps,2, function(x) 10^-(x))
	snps4 <- as.matrix(snps3)
}

snps4 <- apply(snps4,2, function(x) -log10(x))
snps4 <- as.matrix(snps4)
snps4 <- ifelse(snps4 < 0.00001,0.00001,snps4 )
snps4  <-round(snps4, 6)

print("R ALERT: FET p-values combined using Fisher's Method")
outname1 <- paste0(outdir,"/", alertname2, "2.fet")
write.table(snps4, file=outname1, row.names=FALSE, col.names=FALSE)


