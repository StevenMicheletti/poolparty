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

	#turn off sci notation
	options(scipen=999)

print("R ALERT: R structure analysis has begun")

	alertname <- sub('.*/', '', name )
	alertname2 <- sub('\\..*', '', alertname)

#Minor error checking
	if (win.size < 1 | !is.numeric(win.size)) stop ("win.size has to be a positive integer")

#Read freq file produced by sync to af . R, skip the heading 
	infile<-fread(file=name, showProgress=FALSE, header=FALSE)
	infile <- infile[complete.cases(infile), ]
	NPOPS <- ncol(infile) -3
	alert1<- paste0("R ALERT: ", NPOPS, " populations included for structure analyses")
	print(alert1)
	infile$R = rowMaxs(infile[,4:ncol(infile)]) - rowMins(infile[,4:ncol(infile)])

	#total size
	ts <- (max(infile$V2)) + 1000000

#Filter High Allele Diffs
	infileaf <- infile[ which(infile$R > af.filt),]
	infile <- infile[ which(infile$R <= af.filt),]
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

if (win.size > 1) {
  #Calculate SNP density 
  alert5 = paste0("R ALERT: Calculating SNP density at window size of ", win.size, " bp")
  print(alert5)
  denZ <- infile
  denZ$count <- as.numeric(ave(denZ[[1]], denZ$Group, FUN=length))
  denZ$density <- denZ$count / win.size # get proportion of linkage group
  denZ$avpos <-  ((win.size*ceiling(ave(denZ$V2, denZ$Group, FUN=max)/win.size)))
  denZ = denZ[!duplicated(denZ$Group),]
  denzS = 1:nrow(denZ)
  denZI <- as.data.frame(cbind(denZ$V1, denZ$avpos, denzS, denZ$density))
  colnames(denZI) <- c("CHR", "POS", "SNP", "DENSITY")
  denZI$DENSITY<- as.numeric(as.character(denZI$DENSITY))
  minz <- min(denZ$count)
  maxz <- max(denZ$count)
  info <- paste0("R ALERT: Max density is ", maxz, " SNPs per ", win.size, " bp ", "and Min density is ", minz, " SNPs per ", win.size, " bp.")
  print(info)
  info2<- paste0("R ALERT: With a window size of ", win.size, " there are now ", nrow(denZI), " SNP positions")
  print(info2)
  #Write out file
  write.table(denZI, file =tfile, sep='\t', quote=FALSE, col.names= FALSE, row.names = FALSE)
  gpout = paste0(outdir, "/", alertname2, "_density.txt")
  invisible(file.rename(tfile, gpout))
}

print("R ALERT: Processing window size")
	infile$LG=NULL
	infile$V1=NULL

if (win.size > 1) {
  if (method == "mean") {
	infile$V3=NULL
	infile$R=NULL
	infile$V2=NULL
    for (i in 1:(NPOPS)) {
      infile[[i]] <- ave(infile[[i]], infile$Group, FUN=mean)
    }
    infile = infile[!duplicated(infile$Group),]	
    reference = as.data.frame(1:nrow(infile))
    reference$V1 = "p"
	reference=reference$V1
  }
  
  if (method == "random") {
    infile$R=NULL
    infile$V2=NULL
    infile$RAND <- sample(100, size = nrow(infile), replace = TRUE)
    infile <- infile[with(infile, order(Group, -RAND)), ]
    infile = infile[!duplicated(infile$Group),]
    infile$RAND <- NULL
    reference = infile$V3
    infile$V3=NULL
  }
  
  if (method == "rangemax") {
    infile$V2=NULL
    infile <- infile[with(infile, order(Group, -R)), ]
    infile = infile[!duplicated(infile$Group),]
	infile$R=NULL
    reference=infile$V3
    infile$V3=NULL
  }
  
  if (method == "rangemin") {
        infile$V2=NULL
    infile <- infile[with(infile, order(Group, R)), ]
    infile = infile[!duplicated(infile$Group),]
	infile$R=NULL
    reference = infile$V3
        infile$V3=NULL
  }
  
  if (method == "first") {
    infile$R=NULL
    infile <- infile[with(infile, order(Group, V2)), ]
    infile = infile[!duplicated(infile$Group),]
    infile$V2=NULL
    reference = infile$V3
    infile$V3=NULL
  }


  if (method == "last") {
    infile$R=NULL
    infile <- infile[with(infile, order(Group, -V2)), ]
    infile = infile[!duplicated(infile$Group),]
    infile$V2=NULL
    reference = infile$V3
    infile$V3=NULL
  }
rownames(infile) <- infile$Group
infile$Group = NULL
colnames(infile) <- cnames
}

#remove excess columns

if (win.size == 1) {
infile$V3=NULL
infile$V2=NULL
infile$Group = NULL
infile$R = NULL
colnames(infile) <- cnames
}

if (win.size > 1) {
	paleter=paste0("R ALERT: Creating subset frequency table: ", win.size, " bp win size with combine method= ", method)
	print(paleter)
	outnamesub <- paste0(outdir,"/", alertname2, "_subset.fz")
	outfile= as.data.frame(denZI[,1:2])
	outfile$V3=reference
	in2=round(infile,3)
	outfile=cbind(outfile,in2)
	colnames(outfile)[3] <- "Ref"
	write.table(outfile, file=outnamesub , sep='\t',row.names= FALSE, quote=FALSE, col.names=TRUE)
}


print("R ALERT: Creating PCA...")

forPCA <- infile

pca1 = prcomp(forPCA, scale = TRUE)

#Calculate variance explained
pca.eig= pca1$sdev^2
ax1 <-round((pca.eig[1] / sum(pca.eig)*100), digits=1)
ax2 <-round((pca.eig[2] / sum(pca.eig)*100), digits=1)
ax1 <- paste ("PC1 ","(", ax1,"%",")", sep= "")
ax2 <- paste ("PC2 ","(", ax2,"%",")", sep= "")

#Select PC1 and 2 for graphing
plotting = as.data.frame(pca1$rotation[,1:2])

#Make custom Y offset for text
plottingL = plotting
offS = (max(plotting$PC2) - min(plotting$PC2)) * 0.05
plottingL$PC2 = plottingL$PC2 - offS

#make color funcitons
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
bluefunc <- colorRampPalette(c("red","orange","green"))
ctable= as.data.frame(bluefunc(20))

#create sequence bases on 10 color category
coloRZ=seq(from= -.1, to= 1.1, length.out = 21)
ctable$from= coloRZ[1:nrow(ctable)]
ctable$to= coloRZ[2:(nrow(ctable)+1)]

#Apply colors to PC values by bins
plotting2=plotting
plotting2$to =range01(plotting2$PC1)
plotting2$from =range01(plotting2$PC1)
setDT(ctable)
setDT(plotting2)

if (max(plotting2$PC1) - min(plotting2$PC1) > 0) {

	head(plotting2)
	head(ctable)

	setkey(ctable, from, to)

	ans = foverlaps(plotting2, ctable, type="any")

	clist<-as.vector(ans[[1]])
	plotting$clr=paste0(clist,50)
	plotting = as.data.frame(plotting)
	##now Y
	plotting2=plotting
	plotting2$to = range01(plotting2$PC2)
	plotting2$from = range01(plotting2$PC2)
	setDT(ctable)
	setDT(plotting2)
	setkey(ctable, from, to)
	ans = foverlaps(plotting2, ctable, type="any")
	clist<-as.vector(ans[[1]])
	plotting$clr2=paste0(clist,50)
	plotting = as.data.frame(plotting)
	#Expand graph limits by 5%
	expandY = (max(plotting$PC2) - min(plotting$PC2)) * 0.05
	expandX = (max(plotting$PC1) - min(plotting$PC1)) * 0.05
	#Plot 
	outname4 <- paste0(outdir,"/", alertname2, "_PCA.pdf")
	pdf(outname4 , width = 10, height = 8)
	plot(plotting$PC1, plotting$PC2, pch = 20, cex =3, xlab= ax1, ylab =ax2, 
  	   xlim =c(min(plotting$PC1) - expandX, max(plotting$PC1) + expandX ),
    	 ylim =c(min(plotting$PC2) - expandY, max(plotting$PC2) + expandY ),
   	  col=plotting$clr)
  	points(plotting$PC1, plotting$PC2, pch = 20, cex =3,col=plotting$clr2)
  	points(plotting$PC1, plotting$PC2, pch = 20, cex =3,col=plotting$clr2)
  	text(plottingL$PC1,plottingL$PC2, cnames)
	invisible(dev.off())
} else {
	print("R ALERT: No variation, skipping PCA..")
	}

print("R ALERT: Doing the phylogenetics...")

	if (NPOPS <3){
		print("R ALERT: Must have >2 populations to make NJ tree. Skipping..")
	}

if (NPOPS >2){
  print("R ALERT: Proceeding with NJ tree")
affile<-t(infile)

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
  tiplabels(pch=20, col=plotting$clr, cex=3)
  tiplabels(pch=20, col=plotting$cl2, cex=3)
invisible(dev.off())

outname3 <- paste0(outdir,"/", alertname2, "_single.pdf")
 pdf(outname3 , width = 10, height = 8)
  plot(point_est, type = "u")
  nodelabels(bs$BP)
  add.scale.bar()
  tiplabels(pch=20, col=plotting$clr, cex=3)
  tiplabels(pch=20, col=plotting$cl2, cex=3)
invisible(dev.off())
}

print("R ALERT: R NJ analysis done")
