#!bin/R

#Install packages if they don't exist already
	list.of.packages <- c("ggplot2","RColorBrewer","data.table","tidyr")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)) print("RALERT: Installing dependencies for first time use....")
	if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

	

suppressMessages(require(ggplot2))
suppressMessages(require(data.table))
suppressMessages(require(RColorBrewer))
suppressMessages(require(tidyr))

start_time <- Sys.time()

#Get info for bash script
	args <- commandArgs()
	input <- args[6]
	outname <- args[7]
	scaff<- args[8]
	xi <- args[9]
	FUni <- args[10]

if (xi != "NULL") {
	xi <- as.numeric(args[8])
}
	
####LOADING FUNCTIONS from Fariello et al. 2017 ######

# computation of the autocorrelation
autocor=function(x){abs(cor(x[-1],x[-length(x)]))} 

### computation of the lindley process from scores

lindley=function(scores){
  L=length(scores)
	sl=rep(0,L+1)
	for (i in 1:L){
		sl[i+1]= max(0, (sl[i]+scores[i]))
		}
	return(sl[-1])
	}


# computation of the local score from a lindley process and the region under 
# the maximum peak

scorelocal=function(lind){
	nC=max(lind[,1])
	sloc=matrix(0, nrow=nC, ncol=3)
	colnames(sloc)=c('chr','beg','end')
	for (i in 1:nC){
		tmp=which(lind[,1]==i)		
		list=lind[tmp,3]
		M_loc=which.max(list)
		if(length(which(list[1:M_loc]==0))>0){
			m_loc=max(which(list[1:M_loc]==0))+1			
			}else{m_loc=1}
	sloc[i,]=c(i,lind[tmp[m_loc],2],lind[tmp[M_loc],2])	
	}
	return(data.frame(zone=seq(1,nrow(sloc)),sloc))
}

#computation of the significance threshold if the distribution of the p-values is uniform

thresUnif=function(L, cor, xi, alpha = 0.05){
  coefs=list(list('a'=c(-5.5,6.76,-5.66,-2.51),'b'=c(-1.22,3.17,-1.99)),
             list('a'=c(2.47,-4.16,-1.82,-4.58),'b'=c(0.37,2.14,-2.35)),
             list('a'=c(2.04,-5.76,1.04,-6.95),'b'=c(2.55,-0.02,-2.31)),
             list('a'=c(0.22,-4.08,1.16,-9.16),'b'=c(3.45,-0.98,-2.33))
             )
  cors=c(cor^3,cor^2,cor,1)
  if (missing(xi) | !(xi %in% 1:4)){
    print('xi should be 1, 2, 3 or 4') 
    }else{
    a=log(L)+ coefs[[xi]]$a %*%cors
    b=coefs[[xi]]$b %*%cors[-1]
    #then compute the threshold:
  thres = ( log(-log(1-alpha)) - a ) / b
  return(thres)
  }
}

# computation of the significative regions from a lindley process given a significance threshold

sig_sl=function(lind,pos, th){
	zones=c(0,0,0)	
	list=lind
	auxpos=pos
	while(max(list)>=th){
	  M_loc=which.max(list)
		if(length(which(list[1:M_loc]==0))==0){ #the peak is at the beginning of the chrom
			m_loc=1
			zones=rbind(zones, c(auxpos[m_loc],auxpos[M_loc],max(list)))
			tmp=which.min[which(list[M_loc+1:length(list)]==0)] #first 0 score after peak
			list=list[tmp:length(list)]
			auxpos=pos[tmp:length(list)]
			}else{	
				m_loc=max(which(list[1:M_loc]==0))			
				max=max(list)
				zones=rbind(zones, c(auxpos[m_loc+1],auxpos[M_loc],max))
				tmp=which(list[M_loc:length(list)]==0) #first 0 score after peak
				if (length(tmp)>0){
				  auxpos=auxpos[c(1:m_loc,(min(tmp)+M_loc):length(list))]
				  list=list[c(1:m_loc, (min(tmp)+M_loc):length(list))]
				  }else{ #if the peak is at the end of the chromosome
				    auxpos=auxpos[1:m_loc]
				    list=list[1:m_loc]
				    }				
				}
	  }
	zones=matrix(zones, ncol=3)
	zones=data.table(beg=zones[,1],end=zones[,2],peak=zones[,3])
	if (nrow(zones)>1){zones=zones[-1,]}
	return(zones)
	}

### Estimation of Gumble coeficients.

# Estimation of the coeficients of the Gumbel distributions

lineGumb=function(x){
  x1tmp=x
  if (length(table(x1tmp)) >  5){
    Fn = ecdf(x1tmp)
    Fnx=Fn(x1tmp)
    filtre = ( Fnx < 1 - 10**(-1) )  & ( x1tmp > 1 ) #sinon pbs numq
    lm0 = lm(log(-log(Fnx[filtre])) ~ x1tmp[filtre])$coefficients
  }else{lm0=c(0,0)}  
  return(lm0)
}

# Estimation of the coeficients of the polynomes to compute the Gumbel coefficients 
# depending on the length of the chromosomes and the chromosome autocorrelation

coefsGumb=function(mydata, Ls=seq(10000,80000,10000), nSeq=5000, degMod=2){
  bins=seq(0.025,0.975,0.05) 
  coefs=array(0, dim=c(length(bins),5, length(Ls)))  
  coefs[,5,]=matrix(rep(Ls,length(bins)), ncol=length(Ls),byrow=TRUE)
  as=matrix(0, ncol=length(Ls), nrow=length(bins))
  bs=matrix(0, ncol=length(Ls), nrow=length(bins))
  xx=seq(0,max(mydata$lindley), 1)
  for (j in 1:length(Ls)){
    coefs[,1,j]=bins  
    len=Ls[j]
    tmp=sample(nrow(mydata)-len,nSeq,replace=F) #samples the sequences'beginnings
    DTL=data.table(LocScore = vector('numeric', length=nSeq),cors= vector('numeric', length=nSeq)) 
    for (l in 1:nSeq){
      DTL[l]=mydata[seq(tmp[l],length.out=len),.(max(lindley(score)),autocor(pval))]
    }
    DTL[,bin:=which.min(abs(bins-cors)),cors]
    binNE=sort((unique(DTL$bin)))
    coefs[binNE,4,j]=table(DTL$bin)
    DTCoef=DTL[,.(lineGumb(LocScore)),bin]
    coefs[unique(DTL$bin),c(2,3),j]=matrix(DTCoef$V1, ncol=2, byrow=T)
    ys=coefs[,c(2,3),j]%*%rbind(rep(1,length(xx)), xx)
    #pdf(paste('GumbelinesAllxi',xi,'M',len,'.pdf', sep=''))
    #for (i in binNE){
    #  Fn=ecdf(DTL[bin==i,LocScore]) 
    #  tmp=DTL[bin==i,.(LocScore,Fn=Fn(LocScore)),][(Fn< 1 - 10**(-5)) & (LocScore > 1)]
    #  if(i==min(binNE)){
    #    plot(tmp$LocScore, log(-log(tmp$Fn)), col=i,main=paste("M = ",len,sep=""), xlim=quantile(DTL$LocScore,c(0,0.9)),ylim=range(ys[binNE,xx<quantile(DTL$LocScore,0.8)], na.rm=T), pch=20, xlab='Local Score', ylab='log(-log(Fn(LS)))')
    #    }else{points(tmp$LocScore, log(-log(tmp$Fn)), col=i, pch=20)
    #  }
    #  if (sum(ys[i,]!=0)){lines(xx, ys[i,], col=i, lwd=1.5)}
    #}
    #legend("topright",legend=bins[binNE], pch=16,col=seq(1:length(binNE)), title='rho')
    #dev.off()
    as[,j]=as.numeric(coefs[,2,j])
    bs[,j]=as.numeric(coefs[,3,j])		
  }
  
  auxWhich=which(coefs[,4,]> 150)
  rhos=coefs[,1,][auxWhich]
  auxAs=as[auxWhich]-log(coefs[,5,][auxWhich])
  auxBs=bs[auxWhich]
  
  fitA=lm((auxAs)~I(rhos))  
  fitB=lm((auxBs)~I(rhos))  
  
  pdf('FitAandB.pdf')
  par(mfrow=c(1,2))
  plot(rhos, auxAs, pch=16, col=coefs[,5,][auxWhich]/10000, xlab='rho', ylab='a-log(M)')
  xslm=seq(0,1,0.01)
  lines(xslm, fitA$coefficients%*%rbind(rep(1,length(xslm)), xslm))
  plot(rhos, auxBs, pch=16, col=coefs[,5,][auxWhich]/10000, xlab='rho', ylab='b')
  lines(xslm, fitB$coefficients%*%rbind(rep(1,length(xslm)), xslm))
  legend(min(rhos), max(auxBs), legend=sort(unique(coefs[,5,][auxWhich])), col=unique(coefs[,5,][auxWhich]/10000), pch=16, title='M values')
  
  dev.off()
  
  return(list(aCoef=fitA$coefficients, bCoef=fitB$coefficients))
  }

# Computation of the significance threshold given the computed polynomes for computing
# the Gumbel coefficients depending of length and autocorrelation

threshold=function(L, cor, aCoef, bCoef, alpha = 0.05){
  degA=length(aCoef)
  degB=length(bCoef)
  a=log(L)
  b=0
  for (i in 1:degA){
    a=a+aCoef[i]*cor^(i-1)  
    }
  for (i in 1:degB){
    b=b+bCoef[i]*cor^(i-1)  
  }
  #then compute the threshold:
  thres = ( log(-log(1-alpha)) - a ) / b
  return(thres)
  }


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


###### FUNCTIONS DONE #########


alert0=paste0("R ALERT: Loading input data from ", input)
print (alert0)
mydata = fread(input)
colnames(mydata) <- c("chr", "pos", "snp", "pval")

#Remove scaffold headings
mydata <- mydata[!grepl(scaff, mydata$chr),]

#Make chromosomes numbers
#mydata$chr <- as.integer(gsub('[a-zA-Z]', '', mydata$chr ))

Rstore= nrow(mydata)
mydata<- mydata[complete.cases(mydata),]
Xstore= nrow(mydata)

if (Rstore != Xstore) {
	nalert =paste0("R Alert: " ,Rstore-Xstore, " SNPs lost due to NAs in their p-value")
}

if (max(mydata$pval > 1)) {
	print ("R ALERT: P-values are -log 10, converting back to pval")
	mydata$pval <- (10^(-mydata$pval))
}

setkey(mydata, chr)
mydata$pval[mydata$pval==0]=1e-16
Nchr=length(mydata[,unique(chr)])

#Genomewide plots
#This is useful for doing genomewide plots. 

chrInfo=mydata[,.(L=.N,cor=autocor(pval)),chr]
setkey(chrInfo,chr)
data.table(chr=mydata[,unique(chr),], S=cumsum(c(0,chrInfo$L[-Nchr]))) %>% 
  setkey(.,chr) %>% mydata[.,posT:=pos+S]

# To choose the appropriate] threshold ($\xi = 1$ or 2) we look at the distribution of -log10(p value)
print ("R ALERT: Producing p-value distribution and -log10 histogram")
	pname= paste0(outname,"_pvaldist.pdf")
		pdf(pname)
			qplot(pval, data=mydata, geom='histogram', binwidth=0.1, main='P-values histogram')
		invisible(dev.off())

	bname= paste0(outname,"_blog10.pdf")
		pdf(bname)
			qplot(-log10(pval), data=mydata, geom='histogram', binwidth=0.1, main='-log10(p-values) histogram')
		invisible(dev.off())

meanP= mean(-log10(mydata$pval))
stdP= sd(-log10(mydata$pval))

print ("R ALERT: Testing for p-value uniformity with Kolmogorv Smirnov test")
#Kolmogorv Smirnov Test for uniformity

minP=min(mydata$pval)
maxP=max(mydata$pval)

KST=ks.test(unique(mydata$pval), "punif",minP,maxP)[[2]]

if (KST <= 0.05 ) {
	print("R ALERT: Data (p-values) are not uniformly distributed, will apply correction")
}

meanl10<- mean(-log10(mydata$pval))
alert1<- paste0("R ALERT: Mean of log10 is ", round(meanl10,3))
print(alert1)

if (xi == "NULL") {
	print("R ALERT: Determining tuning parameter (xi) using quantiles and mean score")
	# We should verify that the mean of the score is negative. It should be according to the chosen value of xi.
	xi=round((quantile(-log10(mydata$pval), probs = seq(0, 1, 0.05))[[18]]),3)
	alert1<-paste0("R ALERT: Initial tuning parameter set to 85% quantile: ", xi)
	print(alert1)

	mdt=mydata
	mdt[,score:= -log10(pval)-xi]
	meanscore <- mean(mdt$score)
	alert2<-paste0("R ALERT: Initial mean score is: ", meanscore)
	print(alert2)

	if (meanscore >= 0 | meanl10 > 4 ) {
		print ("R ALERT: Mean score is not negative or differentiation is too high trying 90 % quantile")
		mdt2=mydata
		xi=round((quantile(-log10(mydata$pval), probs = seq(0, 1, 0.05))[[19]]),3)
		mdt2[,score:= -log10(pval)-xi]
		meanscore <- mean(mdt2$score)
		alert2<-paste0("R ALERT: Initial mean score is: ", meanscore)
		print(alert2)
			if (meanscore >= 0 | meanl10 > 8) {
					mdt3=mydata
					print ("R ALERT: Mean score is still not negative or differentiation is still too high, trying 95% quantile")
					xi=round((quantile(-log10(mydata$pval), probs = seq(0, 1, 0.05))[[20]]),3)
					mdt3[,score:= -log10(pval)-xi]
					meanscore <- mean(mdt3$score)
					alert2<-paste0("R ALERT: Initial mean score is: ", meanscore)
					print(alert2)
					if (meanscore >= 0) {
						print ("R ALERT: WARNING mean score is not negative. You can specify a custom Xi")
					}
			}
	}
}


mydata[,score:= -log10(pval)-xi]

meanscore <- mean(mydata$score)

alert15<-paste0("R ALERT: Using xi parameter of: ", xi)
alert2<-paste0("R ALERT: Final mean score is: ", meanscore)
print(alert15)
print(alert2)
mydata[,lindley:=lindley(score),chr]

#For uniform distribution

if (KST >= 0.05 | FUni!="NULL" ) {
	xi = round(xi, digits = 0)
	 if (xi > 4){
		xi=4
	 }
	 if (xi < 1) {
		xi=1
	 }
	
	xialtert=paste0("R ALERT: Uniform p-val xi is ", xi)
	print(xialtert)

	print("R ALERT: Calculating significance for a uniform p-val distribution")
	chrInfo[,thG01:=thresUnif(L, cor, xi, 0.01),chr]
	chrInfo[,thG001:=thresUnif(L, cor, xi, 0.001),chr]
	chrInfo[,thG05:=thresUnif(L, cor, xi, 0.05),chr]
	chrInfo[,thG075:=thresUnif(L, cor, xi, 0.075),chr]
	mydata=mydata[chrInfo]

	sigZones05=mydata[,sig_sl(lindley, pos, unique(thG05)),chr]
	sigZones01=mydata[,sig_sl(lindley, pos, unique(thG01)),chr]
	sigZones075=mydata[,sig_sl(lindley, pos, unique(thG075)),chr]
	sigZones001=mydata[,sig_sl(lindley, pos, unique(thG001)),chr]
	
	print("R ALERT: Determining significance threshold and Writing files")
	sig05name = paste0(outname, "_05sig.txt")
	sig01name = paste0(outname, "_01sig.txt")
	sig001name = paste0(outname, "_001sig.txt")
	sig075name = paste0(outname, "_075sig.txt")
	

	ind=which(sigZones05[,peak]>0)
	if (nrow(sigZones05) >0) {
		write.table(sigZones05[ind,],file=sig05name,col.names=T,row.names=F,quote=F)
	}
	ind=which(sigZones01[,peak]>0)
	if (nrow(sigZones01) >0) {
		write.table(sigZones01[ind,],file=sig01name,col.names=T,row.names=F,quote=F)
	}
	ind=which(sigZones001[,peak]>0)
	if (nrow(sigZones001) >0 ) {
		write.table(sigZones001[ind,],file=sig001name,col.names=T,row.names=F,quote=F)
	}
	ind=which(sigZones075[,peak]>0)
	if (nrow(sigZones075) >0 ) {
		write.table(sigZones075[ind,],file=sig075name,col.names=T,row.names=F,quote=F)
	}

	avname = paste0(outname, "_mean_sig.txt")

	#Average of significance thresholds 
	cal05 <- mean(mydata$thG05)
	cal01 <- mean(mydata$thG01)
	cal075 <- mean(mydata$thG075)
	cal001 <- mean(mydata$thG001)

	AVsig=rbind(c(0.001,0.01,0.05,0.075), c(cal001,cal01,cal05,cal075))
	write.table(AVsig,file=avname,col.names=F,row.names=F,quote=F)

	#Chromosome specific thresholds 
	mydatachr <- mydata[,unique(mydata$chr)]
	mydatachr<- mydata[!duplicated(mydata$chr), ]

		cal05chr <- mydatachr$thG05
		cal01chr <- mydatachr$thG01
		cal075chr <- mydatachr$thG075
		cal001chr <- mydatachr$thG001

	chrname = paste0(outname, "_chr_sig.txt")
	CHsig= as.data.frame(cbind(mydatachr$chr,cal001chr,cal01chr,cal05chr,cal075chr))
	colnames(CHsig) <- c("Chr", "sig.001", "sig.01", "sig.05", "sig.075")
	write.table(CHsig,file=chrname,col.names=T,row.names=F,quote=F)


	lsname = paste0(outname, ".ls")
	lsOUT<- as.data.frame(cbind(mydata$chr, mydata$pos, mydata$snp, mydata$lindley))
	write.table(lsOUT,file=lsname,col.names=F,row.names=F,quote=F)
	}else{
	
		print("R ALERT: p-values are not uniform, determining correction coefficients")
		coefsG=coefsGumb(mydata, Ls=seq(30000,60000,10000), nSeq=5000)
		
		print("R ALERT: Determining significance thresholds and writing files")
		chrInfo[,thG05:=threshold(L, cor, coefsG$aCoef, coefsG$bCoef,0.05),]
		chrInfo[,thG01:=threshold(L, cor, coefsG$aCoef, coefsG$bCoef,0.01),]
		chrInfo[,thG001:=threshold(L, cor, coefsG$aCoef, coefsG$bCoef,0.001),]
		chrInfo[,thG075:=threshold(L, cor, coefsG$aCoef, coefsG$bCoef,0.075),]

		mydata=mydata[chrInfo]

		sigZones05=mydata[,sig_sl(lindley, pos, unique(thG05)),chr]
		sigZones01=mydata[,sig_sl(lindley, pos, unique(thG01)),chr]
		sigZones075=mydata[,sig_sl(lindley, pos, unique(thG075)),chr]
		sigZones001=mydata[,sig_sl(lindley, pos, unique(thG001)),chr]

		sig05name = paste0(outname, "_05sig.txt")
		sig01name = paste0(outname, "_01sig.txt")
		sig001name = paste0(outname, "_001sig.txt")
		sig075name = paste0(outname, "_075sig.txt")
	
		ind=which(sigZones05[,peak]>0)
		if (nrow(sigZones05) >0 ) {
			write.table(sigZones05[ind,],file=sig05name,col.names=T,row.names=F,quote=F)
		}
		ind=which(sigZones01[,peak]>0)
		if (nrow(sigZones01) >0 ) {
			write.table(sigZones01[ind,],file=sig01name,col.names=T,row.names=F,quote=F)
		}
		ind=which(sigZones001[,peak]>0)
		if (nrow(sigZones001) >0 ) {
			write.table(sigZones001[ind,],file=sig001name,col.names=T,row.names=F,quote=F)
		}
		ind=which(sigZones075[,peak]>0)
		if (nrow(sigZones075) >0 ) {
			write.table(sigZones075[ind,],file=sig075name,col.names=T,row.names=F,quote=F)
		}
	
		avname = paste0(outname, "_mean_sig.txt")

			#Average of significance thresholds 
			cal05 <- mean(mydata$thG05)
			cal01 <- mean(mydata$thG01)
			cal075 <- mean(mydata$thG075)
			cal001 <- mean(mydata$thG001)
		
		AVsig=rbind(c(0.001,0.01,0.05,0.075), c(cal001,cal01,cal05,cal075))
		write.table(AVsig,file=avname,col.names=F,row.names=F,quote=F)

		#Chromosome specific thresholds 
		mydatachr <- mydata[,unique(mydata$chr)]
		mydatachr<- mydata[!duplicated(mydata$chr), ]

		cal05chr <- mydatachr$thG05
		cal01chr <- mydatachr$thG01
		cal075chr <- mydatachr$thG075
		cal001chr <- mydatachr$thG001


		chrname = paste0(outname, "_chr_sig.txt")
		CHsig<-as.data.frame(cbind(mydatachr$chr,cal001chr,cal01chr,cal05chr,cal075chr))
		colnames(CHsig) <- c("Chr", "sig.001", "sig.01", "sig.05", "sig.075")
		write.table(CHsig,file=chrname,col.names=T,row.names=F,quote=F)

		lsname = paste0(outname, ".ls")
		lsOUT<- as.data.frame(cbind(mydata$chr, mydata$pos, mydata$snp, mydata$lindley))
		write.table(lsOUT,file=lsname,col.names=F,row.names=F,quote=F)

}

end_time <- Sys.time()

timerun = round(difftime(end_time , start_time, units = "mins"),2)
finalaler= paste0("R ALERT: Local.score calculated in ", timerun, " minutes")
print(finalaler)
