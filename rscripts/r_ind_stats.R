#!bin/R

suppressMessages(require('data.table'))
suppressMessages(require('matrixStats'))

#Get info for bash script
args <- commandArgs()
name <- args[6]
outdir <- args[7]

#
outname<- gsub( "RindSTATs.rin", "snp_stats.txt", name )
outname2 <- gsub( "RindSTATs.rin", "ind_stats.txt", name )
alertname <- gsub( "_RindSTATs.rin", "", name )
alertname <- sub('.*/', '', alertname )

#Read in file
filez =fread(file=name, showProgress=FALSE)

#Get populations columns only
pops = as.matrix(filez[,3: ncol(filez)])
#Get population dt with 0s removed
pop0 = pops
is.na(pop0)<- pop0==0

alert1<- paste0("R ALERT: Calculating coverage stats for ", alertname )
print(alert1)

#Get genomic position stats
  #Mean coverage at position [1]
  cov.mat= as.matrix(rowMeans(pops), na.rm=TRUE)
  #Total coverage at position [2]
  cov.mat <- cbind (cov.mat, rowSums(pops))
  #Number of individuals representing position [3]
  cov.mat <- cbind (cov.mat, apply(pops, 1, function(i) sum(i > 0)))
  #Maximum contribution by an individual at position [4]
  cov.mat <- cbind (cov.mat, apply(pops, 1, function(i) max(i)))
  #Maximum % of contribution by an individual [5]
  cov.mat <- cbind (cov.mat, cov.mat[,4] /cov.mat[,2])
  #% of individuals contributing to position [6]
  cov.mat <- cbind (cov.mat, cov.mat[,3] /ncol(pops))
  #Std at position [7]
  cov.mat <- cbind (cov.mat, rowSds(pop0, na.rm=TRUE))
  #Round to reduce file size
  cov.mat <- round(cov.mat, 3)
  
  cov.file <- cbind(filez[,1:2], cov.mat)
  colnames(cov.file) = c("##CHR", "pos","Mean.coverage","Sum.coverage",
                        "Number.individuals", "Maximum.ind", "Maximum.%.ind",
                        "Ind.%", "Std")
                        

alert2<- paste0("R ALERT: Calculating individual stats for ", alertname )
print(alert2)

#Get individual contribution stats 

  #Mean coverage of individual
  ind.cov <- as.data.frame(colMeans(pops))
  #Mean, non-zero coverage of individual
  ind.cov[[2]] <- colMeans(pop0, na.rm=TRUE)
  #Max coverage of this individual
  ind.cov[[3]] <- apply(pops, 2, function(i) max(i))
  #Non-zero standard deviation
  ind.cov[[4]] <- colSds(pop0, na.rm =TRUE)
  #Number of genomic positions with coverage
  ind.cov[[5]] <- apply(pops, 2, function(i) sum(i > 0))
  #Percentage of genome with coverage
  ind.cov[[6]]<- ind.cov[[5]] /nrow(pops)

  ind.cov <- round(ind.cov, 3)
  colnames(ind.cov) = c("##Mean.cov", "Mean.non.zero","Max.cov","std.non.zero",
                         "n.pos", "proportion.SNPs")


##Write out
	alert3<- paste0("R ALERT: Writing output files for ", alertname )
	print(alert3)
	write.table(cov.file, file=outname, sep =" ", row.names = FALSE, col.names = TRUE, quote = FALSE)
	write.table(ind.cov, file=outname2, sep =" ", row.names = FALSE, col.names = TRUE, quote = FALSE)

#DONE