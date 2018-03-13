#!bin/R

print("R ALERT: Checking for R dependencies. Will attempt to install automatically if not present.")

	#Install packages if they don't exist already
	list.of.packages <- c("data.table", "matrixStats","plyr","stringr","ape","fBasics","metap")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)) print("R ALERT: Installing dependencies for first time use....")
	if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

suppressMessages(require(data.table))
suppressMessages(require(matrixStats))
suppressMessages(require(fBasics))
suppressMessages(require(data.table))
suppressMessages(require(ape))
suppressMessages(require(stringr))
suppressMessages(require(plyr))
suppressMessages(require(metap))

args <- commandArgs()
outdir <- args[6]

packageLoaded <- function(name) 0 != length(grep(paste("^package:", 
	name, "$", sep=""), search()))

dtcheck <- packageLoaded("data.table")
mscheck <- packageLoaded("matrixStats")
pcheck <- packageLoaded("plyr")
stcheck <- packageLoaded("stringr")
acheck <- packageLoaded("ape")
fcheck <- packageLoaded("fBasics")
mcheck <- packageLoaded("metap")

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

if (pcheck == FALSE) {
	print("ERROR: R package plyr required! Ensure that this package is installed")
	outname=paste0(outdir,"/","R_ERROR.txt")
	write.table(pcheck, file=outname)
}

if (stcheck == FALSE) {
	print("ERROR: R package stringr required! Ensure that this package is installed")
	outname=paste0(outdir,"/","R_ERROR.txt")
	write.table(stcheck, file=outname)
}

if (acheck == FALSE) {
	print("ERROR: R package ape required! Ensure that this package is installed")
	outname=paste0(outdir,"/","R_ERROR.txt")
	write.table(acheck, file=outname)
}

if (fcheck == FALSE) {
	print("ERROR: R package fBasics required! Ensure that this package is installed")
	outname=paste0(outdir,"/","R_ERROR.txt")
	write.table(fcheck, file=outname)
}

if (mcheck == FALSE) {
	print("ERROR: R package metap required! Ensure that this package is installed")
	outname=paste0(outdir,"/","R_ERROR.txt")
	write.table(mcheck, file=outname)
}
