print("R ALERT: Checking for R dependencies. Will attempt to install automatically if not present.")

	#Install packages if they don't exist already
	list.of.packages <- c("reshape", "fBasics","ggplot2","RColorBrewer")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)) print("R ALERT: Installing dependencies for first time use....")
	if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

 suppressMessages(require(reshape))
 suppressMessages(require(fBasics))
 suppressMessages(require(ggplot2))
 suppressMessages(require(RColorBrewer))


args <- commandArgs()
outdir <- args[6]

packageLoaded <- function(name) 0 != length(grep(paste("^package:", 
	name, "$", sep=""), search()))


dtcheck <- packageLoaded("reshape")
mscheck <- packageLoaded("fBasics")
tcheck <- packageLoaded("ggplot2")
stcheck <- packageLoaded("RColorBrewer")

if (dtcheck  == FALSE) {
	print("ERROR: R package reshape required! Ensure that this package is installed")
	outname=paste0(outdir,"/","R_ERROR.txt")
	write.table(dtcheck, file=outname)
}


if (mscheck == FALSE) {
	print("ERROR: R package fBasics required! Ensure that this package is installed")
	outname=paste0(outdir,"/","R_ERROR.txt")
	write.table(mscheck, file=outname)
}

if (tcheck == FALSE) {
	print("ERROR: R package ggplot2 required! Ensure that this package is installed")
	outname=paste0(outdir,"/","R_ERROR.txt")
	write.table(tcheck, file=outname)
}

if (stcheck == FALSE) {
	print("ERROR: R package RColorBrewer required! Ensure that this package is installed")
	outname=paste0(outdir,"/","R_ERROR.txt")
	write.table(stcheck, file=outname)
}