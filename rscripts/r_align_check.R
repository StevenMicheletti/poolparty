#!bin/R

print("R ALERT: Checking for R dependencies. Will attempt to install automatically if not present.")

	#Install packages if they don't exist already
	list.of.packages <- c("data.table", "matrixStats","tidyr","stringr")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)) print("R ALERT: Installing dependencies for first time use....")
	if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

suppressMessages(require("matrixStats"))
suppressMessages(require("tidyr"))
suppressMessages(require("stringr"))
suppressMessages(require("data.table"))


args <- commandArgs()
outdir <- args[6]


packageLoaded <- function(name) 0 != length(grep(paste("^package:", 
	name, "$", sep=""), search()))


dtcheck <- packageLoaded("data.table")
mscheck <- packageLoaded("matrixStats")
tcheck <- packageLoaded("tidyr")
stcheck <- packageLoaded("stringr")

if (dtcheck  == FALSE) {
	print("ERROR: R package data.table required! Ensure that this package is installed")
	outname=paste0(outdir,"/","R_ERROR.txt")
	write.table(dtcheck, file=outname)
}


if (mscheck == FALSE) {
	print("ERROR: R package matrixStats required! Ensure that this package is installed")
	outname=paste0(outdir,"/","R_ERROR.txt")
	write.table(mscheck, file=outname)
}

if (tcheck == FALSE) {
	print("ERROR: R package tidyr required! Ensure that this package is installed")
	outname=paste0(outdir,"/","R_ERROR.txt")
	write.table(tcheck, file=outname)
}

if (stcheck == FALSE) {
	print("ERROR: R package stringr required! Ensure that this package is installed")
	outname=paste0(outdir,"/","R_ERROR.txt")
	write.table(stcheck, file=outname)
}