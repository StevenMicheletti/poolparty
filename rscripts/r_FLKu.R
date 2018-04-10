#!bin/R

start_time <- Sys.time()

#Install normal packages if they don't exist already
	list.of.packages <- c("ape","phangorn","stringr","data.table")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)) print("R ALERT: Installing dependencies for first time use....")
	if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

#Install biocLite package if it doesn't exist
	list.of.packages <- "qvalue"
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)) print("R ALERT: Installing qvalue from biocLite for first time use....")
	if(length(new.packages)) source ("http://bioconductor.org/biocLite.R") 
	if(length(new.packages)) biocLite("qvalue", suppressUpdates=TRUE)

suppressMessages(require(ape))
suppressMessages(require(qvalue))
suppressMessages(require(data.table))
suppressMessages(require(phangorn))
suppressMessages(library(stringr))


#Get info for bash script
	args <- commandArgs()
	input <- args[6]
	output <- args[7]
	OGr <-args[8]
	
if (OGr == "NULL") {
	OGr=FALSE
}
# Load in FLK.R functions  (https://qgsp.jouy.inra.fr/archives/FLK/FLK.R)

### Compute Reynolds distances for SNP data ####
reynolds=function(Data) 
{
  popnames=rownames(Data)
  npop=nrow(Data) 
  nloc=ncol(Data) 
  dist=matrix(0,nrow=npop,ncol=npop) 
  for (i in 1:(npop-1)) 
    { 
      for (j in (i+1):npop) 
        { 
          pi=Data[i,] 
          pj=Data[j,] 
          dist[i,j]=sum((pi-pj)**2)/sum(pi+pj-2*pi*pj) 
        } 
    } 
  dist=dist+t(dist) 
  rownames(dist)=colnames(dist)=popnames
  return(dist) 
}

midpoint <- function(tree){
  require(phangorn)
    dm = cophenetic(tree)
    tree = unroot(tree)
    rn = max(tree$edge)+1
    maxdm = max(dm)
    ind =  which(dm==maxdm,arr=TRUE)[1,]
    tmproot = Ancestors(tree, ind[1], "parent")
    tree = phangorn:::reroot(tree, tmproot)
    edge = tree$edge
    el = tree$edge.length
    children = tree$edge[,2]
    left = match(ind[1], children)
    tmp = Ancestors(tree, ind[2], "all")
    tmp= c(ind[2], tmp[-length(tmp)])
    right = match(tmp, children)
    if(el[left]>= (maxdm/2)){
         edge = rbind(edge, c(rn, ind[1]))
         edge[left,2] = rn
         el[left] = el[left] - (maxdm/2)
         el = c(el, maxdm/2)
    }
    else{
        sel = cumsum(el[right])
        i = which(sel>(maxdm/2))[1]
        edge = rbind(edge, c(rn, tmp[i]))
        edge[right[i],2] = rn
        eltmp =  sel[i] - (maxdm/2)
#        el = c(el, sel[i] - (maxdm/2))
        el = c(el, el[right[i]] - eltmp)
        el[right[i]] = eltmp
    }
    tree$edge.length = el
    tree$edge=edge
    tree$Nnode  = tree$Nnode+1
    return(phangorn:::reorderPruning(phangorn:::reroot(tree, rn)))
}


### Compute Fij matrix either from SNP Data or from a provided distance
### matrix D, using outgroup as an outgroup population
Fij=function(Data,outgroup=FALSE,D=reynolds(Data))
  {
    require(ape)
    npop=dim(Data)[1]
    D=as.matrix(D)
    PopTree=nj(D)
    if (outgroup==FALSE) {
      PopTree=midpoint(PopTree)
    } else {
      PopTree=root(PopTree,outgroup)
      PopTree=drop.tip(PopTree,outgroup)
      npop=npop-1
    }
    for (i in (which(PopTree$edge.length<0))) {
      PopTree$edge.length[i]=0
    }
    ## get the tree structure as : [ father, son, length]
    edges=cbind(PopTree$edge,PopTree$edge.length)
    #Identify ancestral node as father node with no father itself
    father.nodes=unique(edges[,1])
    son.nodes=unique(edges[,2])
    ancestral=father.nodes[which(is.na(match(father.nodes,son.nodes)))]
    ## Now compute F matrix
    Fij=matrix(0,nrow=npop,ncol=npop,
      dimnames=list(PopTree$tip.label,PopTree$tip.label))
    branch.length=dist.nodes(PopTree)
    route=vector("list",npop)
    tips=1:(npop)
    for (ipop in tips) {
      ## First Fii = distance to the root
      Fij[ipop,ipop]=branch.length[ipop,ancestral]
      ## build route to root
      father=ipop
      while (father!=ancestral) {
        route[[ipop]]=c(route[[ipop]],father)
        edj=which(edges[,2]==father)
        father=edges[edj,1]
      }
      route[[ipop]]=c(route[[ipop]],ancestral)
    }
    ## Now compute the Fij i!= j
    for (ipop in 1:(npop-1)) {
      for (jpop in (ipop+1):(npop)) {
        common.route=intersect(route[[ipop]],route[[jpop]])
        if (length(common.route)>1) {
          for (i in 2:length(common.route)) {
            n1=common.route[i-1]
            n2=common.route[i]
            Fij[ipop,jpop] = Fij[ipop,jpop] + branch.length[n1,n2]
          }
        }
        Fij[jpop,ipop]=Fij[ipop,jpop]
      }
    }
    Fij=2*Fij
    write.table(Fij,"fij.txt",quote=FALSE,col.names=FALSE)
    write(ancestral,file='tree.txt',ncolumns=1)
    write(PopTree$tip.label,file='tree.txt',ncolumns=(npop),append=TRUE)
    write(t(edges),file='tree.txt',ncolumns=3,append=TRUE)
    return(Fij)
  }

### Computes F-LK statistic for a SNP with allele frequenct vector p
### and the inverse of the Fij matrix
FLK.snp=function(p,invF) {
  un=matrix(1,nrow=ncol(invF),ncol=1)
  w=invF%*%un / as.vector( t(un)%*%invF%*%un )
  p0hat=as.double(t(w)%*%p)
  if (p0hat == 0 | p0hat == 1){ 
    Q=0
  }
  else {
    Q=t(p-p0hat*un)%*%invF%*%(p-p0hat*un)
    Q=Q/( p0hat*(1-p0hat) )
    Q=Q * ( 1 - 1/as.double( t(un)%*%invF%*%un ) )
  }
  return(Q)
}

### Computes the numerator of the LK statistic for a SNP with allele
### frequency vector p
LK.num=function(p) {
  if (mean(p)==0 | mean(p)==1) {
    return(0)
  } else {
    return(var(p)/(mean(p)*(1-mean(p))))
  }
}

### Computes the F-LK and original LK test on SNP frequency Data
### conditional and the Fij matrix
FLK=function(Data,Fij) {
  invF=solve(Fij)
  ## consider populations in the Fij matrix only (no outgroup)
  popidx=match(rownames(Fij),rownames(Data))
  freq=Data[popidx,]
  npop=nrow(freq)
  ## F LK test
  F.LK=apply(freq,2,FLK.snp,invF)
  F.LK.p.val=1-sapply(F.LK,pchisq,df=npop-1)
  ## Original LK test
  LK=apply(freq,2,LK.num)
  subset=which(LK!=0)
  LK=(npop-1)*LK/mean(LK[subset])
  LK.p.val=1-sapply(LK,pchisq,df=npop-1)
  pbar=apply(freq,2,mean)
  Ht=2*pbar*(1-pbar)
  return(as.data.frame(cbind(Ht,F.LK,F.LK.p.val,LK,LK.p.val)))
}

##


#Read in allele frequency table (PoolParty format)

inalert=paste0("R ALERT: reading file: ", input)
print(inalert)
 
FLKu <- fread(input, header=T)

#Remove incomplete rows (SNP must have allele frequency in every population)
ALLc <- nrow(FLKu)
FLKu <- FLKu[complete.cases(FLKu),]
NAc <- nrow(FLKu)
	alert0 <- paste0("R ALERT: ", ALLc - NAc, " SNPs removed due to incomplete data")
	print(alert0)
	alert1 <- paste0 ("R ALERT: ", NAc," SNPs will be assessed" )
	print(alert1)

#If unanalyzed format, alter 4th column
if (sum(FLKu[[4]]) < 0) {
	print("R ALERT: Reformatting table")
	FLKu[[4]] =NULL
} 

#Force column names
colnames(FLKu)= c("Chr","Pos", "Ref",  paste0(head(1:(ncol(FLKu)),-3)))


#Create SNP names by joining chr and pos
snpN <- paste0(FLKu$Chr,"_",FLKu$Pos)
FLKu <-FLKu[,4:ncol(FLKu)]


#transpose into FLK input
FLKu <-t(FLKu)
colnames(FLKu)<- snpN

#Run FLK
print ("R ALERT: Calculating Reynolds distances")
DR=reynolds(FLKu)

print ("R ALERT: Running FLK")
F=Fij(FLKu,outgroup=OGr, D=DR)

OutTESTS =FLK(FLKu,F)
qval.FLK <- qvalue(OutTESTS$F.LK.p.val)$qvalues
qval.FLK <- as.data.frame(qval.FLK)
OutTESTS =cbind(OutTESTS,qval.FLK)

split.pos <- as.data.frame(str_split_fixed(rownames(OutTESTS), "_", 2))
colnames(split.pos) = c("Chr", "Pos")

#Formattingoutput
OutTESTS$Ht = round(OutTESTS$Ht,2)
OutTESTS$F.LK = round(OutTESTS$F.LK,2)
OutTESTS$F.LK.p.val = formatC(OutTESTS$F.LK.p.val, format = "e", digits = 3)
OutTESTS$LK = round(OutTESTS$LK,2)
OutTESTS$LK.p.val = formatC(OutTESTS$LK.p.val, format = "e", digits = 3)
OutTESTS$qval.FLK = formatC(OutTESTS$qval.FLK, format = "e", digits = 3)
OutTESTS= cbind(split.pos,OutTESTS)

print ("R ALERT: Writing output table")

#Normal FLK format

falkname1 = paste0(output,"_flk_results.txt")
	write.table(OutTESTS,file=falkname1,col.names=TRUE,row.names=FALSE, quote=FALSE)

end_time <- Sys.time()
timerun = round(difftime(end_time , start_time, units = "mins"),2)
finalaler= paste0("R ALERT: FLK finished in ", timerun, " minutes")
print(finalaler)



