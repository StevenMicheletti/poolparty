###INPUTS####
#############

ppnj <- function (infile,
					win.size = 1,
					method = "none",
					af.filt = 0,
					snp.dens = FALSE,
					chr.num,
					tree.meth = 1,
					boot.s = 100,
					cnames = NULL,
					out.format = "png",
					mask.bs = FALSE
) {


#Dependencies 
require(fBasics)
require(data.table)
require(ape)
require(adegenet)
require(stringr)
require(plyr)


###RUNNING######
###############

#Minor error checking
  if (win.size < 1 | !is.numeric(win.size)) stop ("win.size has to be a positive integer")
  if (af.filt < 0 | af.filt > 1 | !is.numeric(win.size)) stop ("af.filt must be between 0 and 1")
  if (tree.meth < 1 | tree.meth > 5 | !is.numeric(win.size)) stop ("tree.meth must be between 0 and 5")
  if (is.numeric(snp.dens)) stop ("snp.dens should be TRUE or FALSE")
  
#Read freq file produced by sync to af . R

temp4<-fread(infile, skip =1)

NPOPS <- ncol(temp4) -1 

alert1<- paste0(NPOPS, " populations present in this analysis")
print(alert1)

#May be faster as a matrix
##temp4 <- as.matrix(temp4)

#Get allele freq diffs, filter out MAF < 0.05

temp4$R <- rowMaxs(temp4[,2:ncol(temp4)]) - rowMins(temp4[,2:ncol(temp4)])


#Split the names into CHR and POS. NAS will be introduced by scaffolds
headz <- str_split_fixed(temp4$V1, "_", 2)
temp4$POS <- as.numeric (headz[,2])  
temp4$CHR <-(headz[,1])  
temp4$CHR <- as.integer(gsub('[a-zA-Z]', '', temp4$CHR ))


#Add Chromosome back, NAs will be introduced by scaffold
#temp4$CHR <- as.numeric(gsub("chr", "", temp4$CHR))

temp4$CHR[temp4$CHR > chr.num] <- chr.num +1

#Get rid of NAs (usually caused by unanchored regions without positions)
temp4 <- temp4[complete.cases(temp4), ]

#total size
ts <- (max(temp4$POS)) + 1000000

#Filter High Allele Diffs
temp4 <- temp4[ which(temp4$R > af.filt),]

#Get a max position number and break up into linkage groups

if (win.size > 1) {
	roundUP <- function(x, m){x + m - x %% m}
	ts <- roundUP(ts,win.size)
	theseq<- seq(1,ts, by= win.size)
	getname<- (ts/ win.size) -1
	thename<-seq(1,getname, by =1)
	print("Breaking data.frame into linkage groups")
	test <- cut(temp4$POS, breaks = theseq, labels = thename, right = FALSE)
	temp4$LG <- as.data.frame(test)
}

if (win.size == 1) {
	temp4$LG <- 1:nrow(temp4)
}


#Create group name which is CHR_LG
temp4$Group <- paste(temp4$CHR, temp4$LG, sep="_")

#make temp file for mult comps
backup <- temp4

activedir = getwd()
tfile= paste0(activedir,"/","tempaF.tem")

#Remove scaffolds
temp4 <- temp4[ which(temp4$CHR < chr.num + 1),] #get rid of scaffolds


#Calulate SNP density 
if (snp.dens == TRUE ) {
	temp5 <- temp4
	temp5$count <- as.numeric(ave(temp5[[1]], temp5$Group, FUN=length))
	temp5$density <- temp5$count / win.size # get proportion of linkage group
	temp5$avpos <-  ceiling(ave(temp5$POS, temp5$Group, FUN=mean))
	temp5 = temp5[!duplicated(temp5$Group),]
	temp6 <- as.data.frame(cbind(temp5$CHR, temp5$avpos, temp5$density))
	colnames(temp6) <- c("CHR", "POS", "DENSITY")
	minz <- min(temp6$DENSITY)
	maxz <- max(temp6$DENSITY)
	info <- paste0(" Max density is ", maxz, " SNPs per ", win.size, " bp ", "and Min density is ", minz, " SNPs per ", win.size, " bp.")
	print(info)
	write.table(temp6, file =tfile, sep='\t', quote=FALSE, col.names= TRUE, row.names = FALSE)
	gpout = c(infile,".density")
	gpout = paste(gpout, collapse="")
	file.rename(tfile, gpout)
	#make circos density format
		temp7 <- as.data.frame(cbind(temp6$CHR, temp6$POS, temp6$POS, temp6$DENSITY))
		temp7$V1 <- sprintf('ot%i', temp7$V1)
		write.table(temp7, file =tfile, sep='\t', quote=FALSE, col.names= FALSE, row.names = FALSE)
		gpout = c(infile,".density.txt")
		gpout = paste(gpout, collapse="")
		file.rename(tfile, gpout)
}


#Average for linkage group
rownames(temp4) <- temp4$V1
temp4$V1 = NULL

if (win.size > 1) {
	if (method == "mean") {
		for (i in 1:(NPOPS)) {
		temp4[[i]] <- ave(temp4[[i]], temp4$Group, FUN=mean)
		}
	temp4 = temp4[!duplicated(temp4$Group),]	
	}

	if (method == "sdmax") {
		library(matrixStats)
		temp5 <- as.matrix(temp4[,1:NPOPS])
		temp5 <- cbind(temp5,rowSds(temp5))
		temp4$SD <- temp5[, NPOPS+1]
		temp4 <- temp4[with(temp4, order(Group, SD)), ]
		temp4 = temp4[!duplicated(temp4$Group),]
		temp4$SD <- NULL
		detach("package:matrixStats", unload=TRUE)
	}

	if (method == "sdmin") {
		library(matrixStats)
		temp5 <- as.matrix(temp4[,1:NPOPS])
		temp5 <- cbind(temp5,rowSds(temp5))
		temp4$SD <- temp5[, NPOPS+1]
		temp4 <- temp4[with(temp4, order(Group, -SD)), ]
		temp4 = temp4[!duplicated(temp4$Group),]
		temp4$SD <- NULL
		detach("package:matrixStats", unload=TRUE)
	}

	if (method == "random") {
		temp4$RAND <- sample(100, size = nrow(temp4), replace = TRUE)
		temp4 <- temp4[with(temp4, order(Group, -RAND)), ]
		temp4 = temp4[!duplicated(temp4$Group),]
		temp4$RAND <- NULL
	}

	if (method == "rangemax") {
		temp4 <- temp4[with(temp4, order(Group, R)), ]
		temp4 = temp4[!duplicated(temp4$Group),]
	}

	if (method == "rangemin") {
		temp4 <- temp4[with(temp4, order(Group, -R)), ]
		temp4 = temp4[!duplicated(temp4$Group),]
	}

	if (method == "first") {
		temp4 <- temp4[with(temp4, order(Group, -POS)), ]
		temp4 = temp4[!duplicated(temp4$Group),]
	}
}

#remove excess columns
temp4$Group = NULL
temp4$LG = NULL
temp4$CHR = NULL
temp4$POS = NULL
temp4$R = NULL


if (!is.null(cnames)) { 
     colnames(temp4) <- cnames
}

if (is.null(cnames)) { 
colnames(temp4) <- c(as.character(1:NPOPS))
}

temp4<-t(temp4)

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

affile = temp4
#1 = neis
estimate_tr <- function(m) nj(smdist(m, method= tree.meth))

##OR no g dist
#estimate_tr <- function(m) nj(dist(m))

point_est <- estimate_tr(affile)

#Boostrap trees (100)
bs <- boot.phylo(point_est, affile, estimate_tr, trees=TRUE, B= boot.s)  

#Get consensus tree
phy <- ape::consensus(bs$trees, p= .5)

#Get boostrap support values for positions
pp <- prop.part(bs$trees) 
ans <- prop.clades(phy, part = pp, rooted = TRUE) 

#original tree


	if (out.format == "png" | out.format == "all"  ) {
		#plot consensus
		t1name = paste0(infile,"consensus_nj.png")
		png(t1name, width = 10, height = 5, units = 'in', res = 550)
		plot(phy, type = "u")
			if (mask.bs == FALSE){
			nodelabels(ans) }
		dev.off()

		t2name = paste0(infile,"single_nj.png")
		png(t2name, width = 10, height = 5, units = 'in', res = 550)
		#plot point_est with node support
		plot(point_est, type = "u")
			if (mask.bs == FALSE){
			nodelabels(bs$BP) }
		dev.off()

	}

	if (out.format == "pdf" | out.format == "all"  ) {
		#plot consensus
		t1name = paste0(infile,"consensus_nj.pdf")
		pdf(t1name, width = 10, height = 5)
		plot(phy, type = "u")
		add.scale.bar()
			if (mask.bs == FALSE){
			nodelabels(ans) }
		dev.off()

		t2name = paste0(infile,"single_nj.pdf")
		pdf(t2name, width = 10, height = 5)
		#plot point_est with node support
		plot(point_est, type = "u")
		add.scale.bar()
			if (mask.bs == FALSE){
			nodelabels(bs$BP) }
		dev.off()
	}
}

