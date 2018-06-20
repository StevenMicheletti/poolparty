#!bin/R



#PoolParty v0.81
#r_plotter


#plotter.sh pipes arguments
#Uses qqman for quick manhattan plotting

#Install packages if they don't exist already
	list.of.packages <- c("data.table","plyr")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)) print("RALERT: Installing dependencies for first time use....")
	if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

suppressMessages(require(data.table))
suppressMessages(require(plyr))


start_time <- Sys.time()

#Get info for bash script
	args <- commandArgs()
	name <- args[6]
	outname <- args[7]
	outdir <- args[8]
	analtype <- args[9]
	logtrans <- args[10]
	ranges <- args[11]
	scaff <-args[12]
	scalez <-args[13]
	chromosome <-args[14]
	chrrange <- args[15]
	color1 <-args[16]
	color2 <-args[17]
	Gline <-args[18]
	makepdf <-args[19]
	plottype <-args[20]
	minval <-as.numeric(args[21])
        ztrans <- args[22]
	

	filenamezp = paste0(outname,".pdf")
	filenamezn = paste0(outname,".png")
	manename= paste0(outname, " Manhattan Plot")

	
if (scalez == "F" | scalez == "FALSE"  ) {
	scalez="FALSE"
} else {
	scalez ="TRUE"
}

if (logtrans != "NULL") {
	scalez="TRUE"
}

if (Gline != "F") {
	Gline=as.numeric(Gline)
}

if (ztrans != "NULL" & ztrans != "FALSE" & ztrans != "F" ) {
	ztrans="TRUE"
}

	
if (chrrange != "NULL") {
	chrrange1=as.numeric(sub('\\,.*', '', chrrange))
	chrrange2=as.numeric(sub('.*\\,', '', chrrange))
}

if (chromosome != "NULL") {
	chromosome= as.numeric(chromosome)
}


#Uses function from qqman
#qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots
#Turner 2014

if (plottype == "point") {
 manhattan <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", 
     "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
     genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
      ...) 
 {
     CHR = BP = P = index = NULL
     if (!(chr %in% names(x))) 
         stop(paste("Column", chr, "not found!"))
     if (!(bp %in% names(x))) 
         stop(paste("Column", bp, "not found!"))
     if (!(p %in% names(x))) 
         stop(paste("Column", p, "not found!"))
     if (!(snp %in% names(x))) 
         warning(paste("No SNP column found. OK unless you're trying to highlight."))
     if (!is.numeric(x[[chr]])) 
         stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
     if (!is.numeric(x[[bp]])) 
         stop(paste(bp, "column should be numeric."))
     if (!is.numeric(x[[p]])) 
         stop(paste(p, "column should be numeric."))
     d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
     if (!is.null(x[[snp]])) 
         d = transform(d, SNP = x[[snp]])
     d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
     d <- d[order(d$CHR, d$BP), ]
     if (logp) {
         d$logp <- -log10(d$P)
     }
     else {
         d$logp <- d$P
     }
     d$pos = NA
     d$index = NA
     ind = 0
     for (i in unique(d$CHR)) {
         ind = ind + 1
         d[d$CHR == i, ]$index = ind
     }
     nchr = length(unique(d$CHR))
     if (nchr == 1) {
         options(scipen = 999)
         d$pos = d$BP/1e+06
         ticks = floor(length(d$pos))/2 + 1
         xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
         labs = ticks
     }
     else {
         lastbase = 0
         ticks = NULL
         for (i in unique(d$index)) {
             if (i == 1) {
                 d[d$index == i, ]$pos = d[d$index == i, ]$BP
             }
             else {
                 lastbase = lastbase + tail(subset(d, index == 
                   i - 1)$BP, 1)
                 d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
                   lastbase
             }
             ticks = c(ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR == 
                 i, ]$pos))/2 + 1)
         }
         xlabel = "Chromosome"
         labs <- unique(d$CHR)
     }
     xmax = ceiling(max(d$pos) * 1.03)
     xmin = floor(max(d$pos) * -0.03)
     def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
         las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
             ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
     dotargs <- list(...)
     do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
         names(dotargs)]))
     if (!is.null(chrlabs)) {
         if (is.character(chrlabs)) {
             if (length(chrlabs) == length(labs)) {
                 labs <- chrlabs
             }
             else {
                 warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
             }
         }
         else {
             warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
         }
     }
     if (nchr == 1) {
         axis(1, ...)
     }
     else {
         axis(1, at = ticks, labels = labs, ...)
     }
     col = rep(col, max(d$CHR))
     if (nchr == 1) {
         with(d, points(pos, logp, pch = 20, col = col[1], ...))
     }
     else {
         icol = 1
         for (i in unique(d$index)) {
             with(d[d$index == unique(d$index)[i], ], points(pos, 
                 logp, col = col[icol], pch = 20, ...))
             icol = icol + 1
         }
     }
     if (suggestiveline) 
         abline(h = suggestiveline, col = "blue")
     if (genomewideline) 
         abline(h = genomewideline, col = "red")
     if (!is.null(highlight)) {
         if (any(!(highlight %in% d$SNP))) 
             warning("You're trying to highlight SNPs that don't exist in your results.")
         d.highlight = d[which(d$SNP %in% highlight), ]
         with(d.highlight, points(pos, logp, col = "green3", pch = 20, 
             ...))
     }
 }
}

if (plottype == "line") {
 manhattan <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", 
     "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
     genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
      ...) 
 {
     CHR = BP = P = index = NULL
     if (!(chr %in% names(x))) 
         stop(paste("Column", chr, "not found!"))
     if (!(bp %in% names(x))) 
         stop(paste("Column", bp, "not found!"))
     if (!(p %in% names(x))) 
         stop(paste("Column", p, "not found!"))
     if (!(snp %in% names(x))) 
         warning(paste("No SNP column found. OK unless you're trying to highlight."))
     if (!is.numeric(x[[chr]])) 
         stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
     if (!is.numeric(x[[bp]])) 
         stop(paste(bp, "column should be numeric."))
     if (!is.numeric(x[[p]])) 
         stop(paste(p, "column should be numeric."))
     d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
     if (!is.null(x[[snp]])) 
         d = transform(d, SNP = x[[snp]])
     d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
     d <- d[order(d$CHR, d$BP), ]
     if (logp) {
         d$logp <- -log10(d$P)
     }
     else {
         d$logp <- d$P
     }
     d$pos = NA
     d$index = NA
     ind = 0
     for (i in unique(d$CHR)) {
         ind = ind + 1
         d[d$CHR == i, ]$index = ind
     }
     nchr = length(unique(d$CHR))
     if (nchr == 1) {
         options(scipen = 999)
         d$pos = d$BP/1e+06
         ticks = floor(length(d$pos))/2 + 1
         xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
         labs = ticks
     }
     else {
         lastbase = 0
         ticks = NULL
         for (i in unique(d$index)) {
             if (i == 1) {
                 d[d$index == i, ]$pos = d[d$index == i, ]$BP
             }
             else {
                 lastbase = lastbase + tail(subset(d, index == 
                   i - 1)$BP, 1)
                 d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
                   lastbase
             }
             ticks = c(ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR == 
                 i, ]$pos))/2 + 1)
         }
         xlabel = "Chromosome"
         labs <- unique(d$CHR)
     }
     xmax = ceiling(max(d$pos) * 1.03)
     xmin = floor(max(d$pos) * -0.03)
     def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
         las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
             ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
     dotargs <- list(...)
     do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
         names(dotargs)]))
     if (!is.null(chrlabs)) {
         if (is.character(chrlabs)) {
             if (length(chrlabs) == length(labs)) {
                 labs <- chrlabs
             }
             else {
                 warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
             }
         }
         else {
             warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
         }
     }
     if (nchr == 1) {
         axis(1, ...)
     }
     else {
         axis(1, at = ticks, labels = labs, ...)
     }
     col = rep(col, max(d$CHR))
     if (nchr == 1) {
         with(d, points(pos, type="l", logp, pch = 20, col = col[1], ...))
     }
     else {
         icol = 1
         for (i in unique(d$index)) {
             with(d[d$index == unique(d$index)[i], ], points(pos, type="l",
                 logp, col = col[icol], pch = 20, ...))
             icol = icol + 1
         }
     }
     if (suggestiveline) 
         abline(h = suggestiveline, col = "blue")
     if (genomewideline) 
         abline(h = genomewideline, col = "red")
     if (!is.null(highlight)) {
         if (any(!(highlight %in% d$SNP))) 
             warning("You're trying to highlight SNPs that don't exist in your results.")
         d.highlight = d[which(d$SNP %in% highlight), ]
         with(d.highlight, points(pos, logp, col = "green3", pch = 20, 
             ...))
     }
 }
}

if (plottype == "bar") {
 manhattan <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", 
     "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
     genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
      ...) 
 {
     CHR = BP = P = index = NULL
     if (!(chr %in% names(x))) 
         stop(paste("Column", chr, "not found!"))
     if (!(bp %in% names(x))) 
         stop(paste("Column", bp, "not found!"))
     if (!(p %in% names(x))) 
         stop(paste("Column", p, "not found!"))
     if (!(snp %in% names(x))) 
         warning(paste("No SNP column found. OK unless you're trying to highlight."))
     if (!is.numeric(x[[chr]])) 
         stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
     if (!is.numeric(x[[bp]])) 
         stop(paste(bp, "column should be numeric."))
     if (!is.numeric(x[[p]])) 
         stop(paste(p, "column should be numeric."))
     d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
     if (!is.null(x[[snp]])) 
         d = transform(d, SNP = x[[snp]])
     d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
     d <- d[order(d$CHR, d$BP), ]
     if (logp) {
         d$logp <- -log10(d$P)
     }
     else {
         d$logp <- d$P
     }
     d$pos = NA
     d$index = NA
     ind = 0
     for (i in unique(d$CHR)) {
         ind = ind + 1
         d[d$CHR == i, ]$index = ind
     }
     nchr = length(unique(d$CHR))
     if (nchr == 1) {
         options(scipen = 999)
         d$pos = d$BP/1e+06
         ticks = floor(length(d$pos))/2 + 1
         xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
         labs = ticks
     }
     else {
         lastbase = 0
         ticks = NULL
         for (i in unique(d$index)) {
             if (i == 1) {
                 d[d$index == i, ]$pos = d[d$index == i, ]$BP
             }
             else {
                 lastbase = lastbase + tail(subset(d, index == 
                   i - 1)$BP, 1)
                 d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
                   lastbase
             }
             ticks = c(ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR == 
                 i, ]$pos))/2 + 1)
         }
         xlabel = "Chromosome"
         labs <- unique(d$CHR)
     }
     xmax = ceiling(max(d$pos) * 1.03)
     xmin = floor(max(d$pos) * -0.03)
     def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
         las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
             ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
     dotargs <- list(...)
     do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
         names(dotargs)]))
     if (!is.null(chrlabs)) {
         if (is.character(chrlabs)) {
             if (length(chrlabs) == length(labs)) {
                 labs <- chrlabs
             }
             else {
                 warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
             }
         }
         else {
             warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
         }
     }
     if (nchr == 1) {
         axis(1, ...)
     }
     else {
         axis(1, at = ticks, labels = labs, ...)
     }
     col = rep(col, max(d$CHR))
     if (nchr == 1) {
         with(d, points(pos, type="h", logp, pch = 20, col = col[1], ...))
     }
     else {
         icol = 1
         for (i in unique(d$index)) {
             with(d[d$index == unique(d$index)[i], ], points(pos, type="h",
                 logp, col = col[icol], pch = 20, ...))
             icol = icol + 1
         }
     }
     if (suggestiveline) 
         abline(h = suggestiveline, col = "blue")
     if (genomewideline) 
         abline(h = genomewideline, col = "red")
     if (!is.null(highlight)) {
         if (any(!(highlight %in% d$SNP))) 
             warning("You're trying to highlight SNPs that don't exist in your results.")
         d.highlight = d[which(d$SNP %in% highlight), ]
         with(d.highlight, points(pos, logp, col = "green3", pch = 20, 
             ...))
     }
 }
}


#If ranges is available, combine CHR file
if (ranges != "NULL") { 
	print("RALERT: Reading files")
    chrm <- read.delim(ranges, as.is=TRUE, header=TRUE)
	colnames(chrm) <- c("CHR","BP")
	chrm$CHR <- as.integer(gsub('[a-zA-Z]', '', chrm$CHR ))
	chrm$BP <- as.integer(gsub('[a-zA-Z]', '', chrm$BP ))
    snps <- fread(name, stringsAsFactors=FALSE)
    snps <- snps [complete.cases(snps), ]
	print("RALERT: Formatting files")
    colnames(snps) <- c("CHR", "BP","SNP","FST")
	snps <- snps[!grepl(scaff, snps $CHR),]
    snps$CHR <- as.integer(gsub('[a-zA-Z]', '', snps$CHR ))
    snps$BP <- as.integer(gsub('[a-zA-Z]', '', snps$BP ))
    snps$FST[snps$FST==0]=1e-16
	snps<- rbind.fill(snps, chrm)
	if (logtrans == "NULL") {
		snps[is.na(snps)] <- 0
	}
	if (logtrans != "NULL") {
		snps[is.na(snps)] <- .9
	}
	snps <- snps[order(snps$CHR, snps$BP),]
	snps$SNP <- 1:nrow(snps) 
	alertz= paste0("RALERT: There are ", nrow(snps), " SNP positions retained for plotting")
	print(alertz)
      if (logtrans != "NULL") {
        snps$FST <- -log10(snps$FST)
	    snps[which(!is.finite(snps$FST))] <- 0
      }
      if (ztrans == "TRUE") {
        snps$FST <- scale(snps$FST)
	    snps[which(!is.finite(snps$FST))] <- 0
      }

	snps <- snps[(snps$FST >=minval),]
	#use scale, use chromosome, use ranges 
	CHRNz <- as.numeric(length(unique(snps$CHR)))
	chrt<- paste0 ("RALERT: Creating plots from genome with ", CHRNz, " chromosomes")
	print(chrt)

	if (scalez == "TRUE" & chromosome !="NULL" & chrrange !="NULL") {
		highY= (max(snps$FST) + max(snps$FST) * 0.10)
		xlabler= paste0( "Position on chromosome ", chromosome, " (mb)")
		if (makepdf != "NULL") {
			pdf(filenamezp, width = 10, height = 5)
				manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(chrrange1,chrrange2), ylim=c(0, highY), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
				main= manename )
		}
		invisible(dev.off())
		png(filenamezn, width = 10, height = 5 , units = 'in', res = 550)
			manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(chrrange1,chrrange2), ylim=c(0, highY), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
			main= manename )
		invisible(dev.off())
	}
	#don't use scale, don't use chromosome, don't use ranges 
	if (scalez == "FALSE"  & chromosome =="NULL" & chrrange =="NULL") {
		if (makepdf != "NULL") {
			pdf(filenamezp, width = 10, height = 5)
				manhattan(snps, suggestiveline = F, genomewideline = Gline,  ylim=c(0, 1.05), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= "Chromosome", 
				main= manename )
			invisible(dev.off())
		}
		png(filenamezn, width = 10, height = 5 , units = 'in', res = 550)
			manhattan(snps, suggestiveline = F, genomewideline = Gline,  ylim=c(0, 1.05), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= "Chromosome", 
			main= manename )
		invisible(dev.off())
	}
	#don't use scale, use chromosome, don't use ranges
	if (scalez == "FALSE" & chromosome !="NULL" & chrrange =="NULL") {
		sVAL <- subset(snps, CHR == chromosome)
		xlabler= paste0( "Position on chromosome ", chromosome, " (mb)")
		if (makepdf != "NULL") {
			pdf(filenamezp, width = 10, height = 5)
				manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(0,max(sVAL$BP)/1000000), ylim=c(0, 1.05),  p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
				main= manename )
			invisible(dev.off())
		}
		png(filenamezn, width = 10, height = 5 , units = 'in', res = 550)
			manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(0,max(sVAL$BP)/1000000), ylim=c(0, 1.05), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
			main= manename )
		invisible(dev.off())
	}
	#Use scale, use chromosome, don't use ranges
	if (scalez == "TRUE" & chromosome !="NULL" & chrrange =="NULL") {
		sVAL <- subset(snps, CHR == chromosome)
		xlabler= paste0( "Position on chromosome ", chromosome, " (mb)")
		highY= (max(snps$FST) + max(snps$FST) * 0.10)
		if (makepdf != "NULL") {
			pdf(filenamezp, width = 10, height = 5)
				manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(0,max(sVAL$BP)/1000000), ylim=c(0, highY), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
				main= manename )
			invisible(dev.off())
		}
		png(filenamezn, width = 10, height = 5 , units = 'in', res = 550)
			manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(0,max(sVAL$BP)/1000000), ylim=c(0, highY), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
			main= manename )
		invisible(dev.off())
	}
	
	#Don't Use scale,  use chromosome, use ranges
	if (scalez == "FALSE" & chromosome !="NULL" & chrrange !="NULL") {
		xlabler= paste0( "Position on chromosome ", chromosome, " (mb)")
		if (makepdf != "NULL") {
			pdf(filenamezp, width = 10, height = 5)
				manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(chrrange1,chrrange2), ylim=c(0, 1.05), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
				main= manename )
			invisible(dev.off())
		}
		png(filenamezn, width = 10, height = 5 , units = 'in', res = 550)
			manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(chrrange1,chrrange2), ylim=c(0, 1.05), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
			main= manename )
		invisible(dev.off())
	}
	
	#Use scale,  don't use chromosome, dont use ranges #COMPLETE
	if (scalez == "TRUE" & chromosome =="NULL" & chrrange =="NULL") {
		highY= (max(snps$FST) + max(snps$FST) * 0.10)
		if (makepdf != "NULL") {
			pdf(filenamezp, width = 10, height = 5)
				manhattan(snps, suggestiveline = F, genomewideline = Gline, ylim=c(0, highY), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= "Chromosome", 
				main= manename )
			invisible(dev.off())
		}
		png(filenamezn, width = 10, height = 5 , units = 'in', res = 550)
			manhattan(snps, suggestiveline = F, genomewideline = Gline, ylim=c(0, highY), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= "Chromosome", 
			main= manename )
		invisible(dev.off())
	}
}

#If ranges is not available, use data extent 
if (ranges == "NULL") { 
	print("RALERT: Reading files")
    snps <- fread(name, stringsAsFactors=FALSE)
    snps <- snps [complete.cases(snps), ]
	print("RALERT: Formatting files")
    colnames(snps) <- c("CHR", "BP","SNP","FST")
	snps <- snps[!grepl(scaff, snps $CHR),]
    snps$CHR <- as.integer(gsub('[a-zA-Z]', '', snps$CHR ))
    snps$BP <- as.integer(gsub('[a-zA-Z]', '', snps$BP ))
    snps$FST[snps$FST==0]=1e-16
	snps[is.na(snps)] <- 0
	snps$SNP <- 1:nrow(snps) 
	alertz= paste0("RALERT: There are ", nrow(snps), " SNP positions retained for plotting")
	print(alertz)
      if (logtrans != "NULL") {
        snps$FST <- -log10(snps$FST)
	snps[which(!is.finite(snps$FST))] <- 0
      }
      if (ztrans == "TRUE") {
        snps$FST <- scale(snps$FST)
	snps[which(!is.finite(snps$FST))] <- 0
      }

	snps <- snps[(snps$FST >= minval),]
	#use scale, use chromosome, use ranges 
	CHRNz <- as.numeric(length(unique(snps$CHR)))
	chrt<- paste0 ("RALERT: Creating plots from genome with ", CHRNz, " chromosomes")
	print(chrt)
	if (scalez == "TRUE" & chromosome !="NULL" & chrrange !="NULL") {
		highY= (max(snps$FST) + max(snps$FST) * 0.10)
		xlabler= paste0( "Position on chromosome ", chromosome, " (mb)")
		if (makepdf != "NULL") {
			pdf(filenamezp, width = 10, height = 5)
				manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(chrrange1,chrrange2), ylim=c(0, highY), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
				main= manename )
			invisible(dev.off())
		}
		png(filenamezn, width = 10, height = 5 , units = 'in', res = 550)
			manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(chrrange1,chrrange2), ylim=c(0, highY), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
			main= manename )
		invisible(dev.off())
	}
	#don't use scale, don't use chromosome, don't use ranges 
	if (scalez == "FALSE" & chromosome =="NULL" & chrrange =="NULL") {
		if (makepdf != "NULL") {
			pdf(filenamezp, width = 10, height = 5)
				manhattan(snps, suggestiveline = F, genomewideline = Gline, ylim=c(0, 1.05), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= "Chromosome", 
				main= manename )
			invisible(dev.off())
		}
		png(filenamezn, width = 10, height = 5 , units = 'in', res = 550)
			manhattan(snps, suggestiveline = F, genomewideline = Gline,  ylim=c(0, 1.05), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= "Chromosome", 
			main= manename )
		invisible(dev.off())
	}
	#don't use scale, use chromosome, don't use ranges
	if (scalez == "FALSE" & chromosome !="NULL" & chrrange =="NULL") {
		xlabler= paste0( "Position on chromosome ", chromosome, " (mb)")
		sVAL <- subset(snps, CHR == chromosome)
		if (makepdf != "NULL") {
			pdf(filenamezp, width = 10, height = 5)
				manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(0,max(sVAL$BP)/1000000), ylim=c(0, 1.05), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
				main= manename )
			invisible(dev.off())
		}
		png(filenamezn, width = 10, height = 5 , units = 'in', res = 550)
			manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(0,max(sVAL$BP)/1000000), ylim=c(0, 1.05), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
			main= manename )
		invisible(dev.off())
	}
	#Use scale, use chromosome, don't use ranges
	if (scalez == "TRUE" & chromosome !="NULL" & chrrange =="NULL") {
		xlabler= paste0( "Position on chromosome ", chromosome, " (mb)")
		highY= (max(snps$FST) + max(snps$FST) * 0.10)
		sVAL <- subset(snps, CHR == chromosome)
		if (makepdf != "NULL") {
			pdf(filenamezp, width = 10, height = 5)
				manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(0,max(sVAL$BP)/1000000), ylim=c(0, highY), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
				main= manename )
			invisible(dev.off())
		}
		png(filenamezn, width = 10, height = 5 , units = 'in', res = 550)
			manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(0,max(sVAL$BP)/1000000), ylim=c(0, highY), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
			main= manename )
		invisible(dev.off())
	}
	
	#Don't Use scale,  use chromosome, use ranges
	if (scalez == "FALSE" & chromosome !="NULL" & chrrange !="NULL") {
		xlabler= paste0( "Position on chromosome ", chromosome, " (mb)")
		if (makepdf != "NULL") {
			pdf(filenamezp, width = 10, height = 5)
				manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(chrrange1,chrrange2), ylim=c(0, 1.05), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
				main= manename )
			invisible(dev.off())
		}
		png(filenamezn, width = 10, height = 5 , units = 'in', res = 550)
			manhattan(subset(snps, CHR==chromosome), suggestiveline = F, genomewideline = Gline, xlim=c(chrrange1,chrrange2), ylim=c(0, 1.05), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= xlabler, 
			main= manename )
		invisible(dev.off())
	}
	
	#Use scale,  don't use chromosome, dont use ranges
	if (scalez == "TRUE" & chromosome =="NULL" & chrrange =="NULL") {
		highY= (max(snps$FST) + max(snps$FST) * 0.10)
		if (makepdf != "NULL") {
			pdf(filenamezp, width = 10, height = 5)
				manhattan(snps, suggestiveline = F, genomewideline = Gline, ylim=c(0, highY), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= "Chromosome", 
				main= manename )
			invisible(dev.off())
		}
		png(filenamezn, width = 10, height = 5 , units = 'in', res = 550)
			manhattan(snps, suggestiveline = F, genomewideline = Gline, ylim=c(0, highY), p="FST", logp=FALSE, col = c(color1, color2), cex=0.1, cex.axis = 0.3, ylab=analtype , xlab= "Chromosome", 
			main= manename )
		invisible(dev.off())
	}
}


end_time <- Sys.time()
timerun = round(difftime(end_time , start_time, units = "secs"),2)

finalaler= paste0("RALERT: Plots created in ", timerun, " seconds")
print(finalaler)
