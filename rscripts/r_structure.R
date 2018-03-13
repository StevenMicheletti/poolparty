#!bin/R
suppressMessages(require(fBasics))
suppressMessages(require(data.table))
suppressMessages(require(ape))
suppressMessages(require(stringr))
suppressMessages(require(plyr))

#pass arguments
	args <- commandArgs()
	name <- args[6]
	outdir <- args[7]
	#AF filter
	af.filt <- as.numeric(args[8])
	#win.size
	win.size <- as.numeric(args[9])
	#method
	method <- args[10]
	#bootstraps
	boot.s <- as.numeric(args[11])
	#cnamefile
	cnamefile <- args[12]

	#Grab column names
	cnames<-as.character(read.table(file=cnamefile , header=F))

print("R ALERT: R structure analysis has begun")

	alertname <- sub('.*/', '', name )
	alertname2 <- sub('\\..*', '', alertname)

	#Minor error checking
	if (win.size < 1 | !is.numeric(win.size)) stop ("win.size has to be a positive integer")

#Read freq file produced by sync to af . R, skip the heading 
	infile<-fread(file=name, showProgress=FALSE)
	infile <- infile[complete.cases(infile), ]
	NPOPS <- ncol(infile) -3
	alert1<- paste0("R ALERT: ", NPOPS, " populations included for structure analyses")
	print(alert1)
	infile$R = rowMaxs(infile[,4:ncol(infile)]) - rowMins(infile[,4:ncol(infile)])

	#total size
	ts <- (max(infile$V2)) + 1000000

#Filter High Allele Diffs
	infileaf <- infile[ which(infile$R > af.filt),]
	infile <- infile[ which(infile$R < af.filt),]
	alert3 <- paste0("R ALERT: ", nrow(infileaf), " SNPs removed due to Maximum allele frequency filter of: ", af.filt)
	alert4 <- paste0("R ALERT: ", nrow(infile), " SNPs still remain")
	print(alert3)
	print(alert4)

#Get a max position number and break up into linkage groups
	if (win.size > 1) {
		roundUP <- function(x, m){x + m - x %% m}
		ts <- roundUP(ts,win.size)
		theseq<- seq(1,ts, by= win.size)
		getname<- (ts/ win.size) -1
		thename<-seq(1,getname, by =1)
		print("R ALERT: Breaking data.frame into windowed linkage groups")
		test <- cut(infile$V2, breaks = theseq, labels = thename, right = FALSE)
		infile$LG <- as.data.frame(test)
	}

	if (win.size == 1) {
		infile$LG <- 1:nrow(infile)
	}

#Create group name which is CHR_LG
	infile$Group <- paste(infile$V1, infile$LG, sep="_")

#Create output name
tfile=paste0(outdir, "/", alertname2, "_TEMP")
#Calculate SNP density 
  alert5 = paste0("R ALERT: Calculating SNP density at window size of ", win.size, " bp")
  print(alert5)
  denZ <- infile
  denZ$count <- as.numeric(ave(denZ[[1]], denZ$Group, FUN=length))
  denZ$density <- denZ$count / win.size # get proportion of linkage group
  denZ$avpos <-  ceiling(ave(denZ$V2, denZ$Group, FUN=mean))
  denZ = denZ[!duplicated(denZ$Group),]
  denZI <- as.data.frame(cbind(denZ$V1, denZ$avpos, denZ$density))
  colnames(denZI) <- c("CHR", "POS", "DENSITY")
  denZI$DENSITY<- as.numeric(as.character(denZI$DENSITY))
  minz <- min(denZI$DENSITY)
  maxz <- max(denZI$DENSITY)
  info <- paste0("R ALERT: Max density is ", maxz, " SNPs per ", win.size, " bp ", "and Min density is ", minz, " SNPs per ", win.size, " bp.")
  print(info)
  #Write out file
  write.table(denZI, file =tfile, sep='\t', quote=FALSE, col.names= TRUE, row.names = FALSE)
  gpout = paste0(outdir, "/", alertname2, "_density.txt")
  invisible(file.rename(tfile, gpout))

#remove excess columns
	infile$V3=NULL
	infile$LG = NULL
	infile$V1 = NULL
	infile$V2 = NULL
	infile$R = NULL

print("R ALERT: Processing window size")

if (win.size > 1) {
  if (method == "mean") {
    for (i in 1:(NPOPS)) {
      infile[[i]] <- ave(infile[[i]], infile$Group, FUN=mean)
    }
    infile = infile[!duplicated(infile$Group),]	
  }
if (method == "sdmax") {
    library(matrixStats)
    temp5 <- as.matrix(infile[,1:NPOPS])
    temp5 <- cbind(temp5,rowSds(temp5))
    infile$SD <- temp5[, NPOPS+1]
    infile <- infile[with(infile, order(Group, SD)), ]
    infile = infile[!duplicated(infile$Group),]
    infile$SD <- NULL
    detach("package:matrixStats", unload=TRUE)
  }
  
  if (method == "sdmin") {
    library(matrixStats)
    temp5 <- as.matrix(infile[,1:NPOPS])
    temp5 <- cbind(temp5,rowSds(temp5))
    infile$SD <- temp5[, NPOPS+1]
    infile <- infile[with(infile, order(Group, -SD)), ]
    infile = infile[!duplicated(infile$Group),]
    infile$SD <- NULL
    detach("package:matrixStats", unload=TRUE)
  }
  
  if (method == "random") {
    infile$RAND <- sample(100, size = nrow(infile), replace = TRUE)
    infile <- infile[with(infile, order(Group, -RAND)), ]
    infile = infile[!duplicated(infile$Group),]
    infile$RAND <- NULL
  }
  
  if (method == "rangemax") {
    infile <- infile[with(infile, order(Group, R)), ]
    infile = infile[!duplicated(infile$Group),]
  }
  
  if (method == "rangemin") {
    infile <- infile[with(infile, order(Group, -R)), ]
    infile = infile[!duplicated(infile$Group),]
  }
  
  if (method == "first") {
    infile <- infile[with(infile, order(Group, -POS)), ]
    infile = infile[!duplicated(infile$Group),]
  }
  if (method == "none") {
    infile <- infile[with(infile, order(Group, -POS)), ]
    infile = infile[!duplicated(infile$Group),]
  }
}


print("R ALERT: Doing the phylogenetics...")
rownames(infile) <- infile$Group
infile$Group = NULL
colnames(infile) <- cnames

	if (NPOPS <3){
		print("R ALERT: Must have >2 populations to make NJ tree. Skipping..")
	}

if (NPOPS >2){
  print("R ALERT: Proceeding with NJ tree")
affile<-t(infile)
options(warn=-1)
# Get distance tree (manhattan = absolute distance between vectors)

smdist<- function (x, method = 1, diag = FALSE, upper = FALSE) 
{
  METHODS = c("Nei", "Edwards", "Reynolds", "Rodgers", "Provesti")
  if (all((1:5) != method)) {
    cat("1 = Nei 1972\n")
    cat("2 = Edwards 1971\n")
    cat("3 = Reynolds, Weir and Coockerman 1983\n")
    cat("4 = Rodgers 1972\n")
    cat("5 = Provesti 1975\n")
    method <- as.integer(readLines(n = 1))
  }
  
  if (all((1:5) != method)) 
    (stop("Non convenient method number"))
  
  nloc <- ncol(x)
  nlig <- nrow(x)
  loc.fac <- nloc
  if (method == 1) {
    d <- x %*% t(x)
    vec <- sqrt(diag(d))
    d <- d/vec[col(d)]
    d <- d/vec[row(d)]
    d <- -log(d)
    d <- as.dist(d)
  }
  else if (method == 2) {
    x <- sqrt(x)
    d <- x %*% t(x)
    d <- 1 - d/nloc
    diag(d) <- 0
    d <- sqrt(d)
    d <- as.dist(d)
  }
  else if (method == 3) {
    denomi <- x %*% t(x)
    vec <- apply(x, 1, function(x) sum(x * x))
    d <- -2 * denomi + vec[col(denomi)] + vec[row(denomi)]
    diag(d) <- 0
    denomi <- 2 * nloc - 2 * denomi
    diag(denomi) <- 1
    d <- d/denomi
    d <- sqrt(d)
    d <- as.dist(d)
  }
  else if (method == 4) {
    kx <- lapply(split(x, loc.fac[col(x)]), matrix, nrow = nlig)
    dcano <- function(mat) {
      daux <- mat %*% t(mat)
      vec <- diag(daux)
      daux <- -2 * daux + vec[col(daux)] + vec[row(daux)]
      diag(daux) <- 0
      daux <- sqrt(0.5 * daux)
      return(daux)
    }
    d <- matrix(0, nlig, nlig)
    for (i in 1:length(kx)) {
      d <- d + dcano(kx[[i]])
    }
    d <- d/length(kx)
    d <- as.dist(d)
  }
  
  attr(d, "Size") <- nlig
  attr(d, "Labels") <- rownames(x)
  attr(d, "Diag") <- diag
  attr(d, "Upper") <- upper
  attr(d, "method") <- METHODS[method]
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  return(d)
}

options(warn=0)

estimate_tr <- function(m) nj(smdist(m, method= 1))

point_est <- estimate_tr(affile)

#Boostrap trees (100)
bs <- boot.phylo(point_est, affile, estimate_tr, trees=TRUE, B= boot.s)  

#Get consensus tree
phy <- ape::consensus(bs$trees, p= .5)
#Get boostrap support values for positions
pp <- prop.part(bs$trees) 
ans <- prop.clades(phy, part = pp, rooted = TRUE) 

#plot consensus
outname2 <- paste0(outdir,"/", alertname2, "_consensus.pdf")
  pdf(outname2, width = 10, height = 8)
  plot(phy, type = "u")
  nodelabels(bs$BP)
  add.scale.bar()
invisible(dev.off())

outname3 <- paste0(outdir,"/", alertname2, "_single.pdf")
  pdf(outname3 , width = 10, height = 8)
  plot(point_est, type = "u")
  nodelabels(bs$BP)
  add.scale.bar()
invisible(dev.off())
}

print("R ALERT: R NJ analysis done")
