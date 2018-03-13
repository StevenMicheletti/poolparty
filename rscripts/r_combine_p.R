#!bin/R
suppressMessages(require(data.table))
suppressMessages(require(metap))

args <- commandArgs()
name <- args[6]
outdir <- args[7]

	alertname <- sub('.*/', '', name )
	alertname2 <- sub('\\..*', '', alertname)

snps <- fread(name, stringsAsFactors=FALSE, header=F, skip=1)
snps <- as.data.frame(snps)
snps <- as.data.frame(snps[,colSums(is.na(snps)) == 0])

print(ncol(snps))

#Get p-val for each column
if (ncol(snps) > 1) {
	snps3 <- apply(snps,2, function(x) 10^-(x))
	snps4 <- as.matrix(apply(snps3,1, function(x) sumlog(x)[[3]])) 
}
if (ncol(snps) < 2) {
	snps3 <- apply(snps,2, function(x) 10^-(x))
	snps4 <- as.matrix(snps3)
}

print("R ALERT: FET p-values combined using Fisher's Method")
outname1 <- paste0(outdir,"/", alertname2, "2.fet")
write.table(snps4, file=outname1, row.names=FALSE, col.names=FALSE)