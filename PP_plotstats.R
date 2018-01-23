#Install packages if they don't exist already
#PP_plotstats.R by Steven Micheletti
#Version B1 1/23/2018


list.of.packages <- c("reshape", "fBasics","ggplot2","RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) print("installing dependencies for first time use....")
if(length(new.packages)) install.packages(new.packages)

plot.stats <- function (infile,
                       outdir=NULL,
                      outfile="Plot_out"
) {
  
  
  #Required Packages 
  require(reshape)
  require(fBasics)
  require(ggplot2)
  require(RColorBrewer)

  infile <- read.delim(infile, header= F, skip=1, sep = " ")
  outname = outfile

  currentD <- getwd()
  if (is.null(outdir)) {
    a1 <- paste0("Files will be written to ", getwd())
    print(a1)
  }
    
  if (!is.null(outdir)) {
    dir.create(file.path(outdir), showWarnings = FALSE)
    setwd(file.path(outdir))
    a1 <- paste0("Files will be written to ", outdir)
    print(a1)
  }
  
  
  infile <- infile[,colSums(is.na(infile))<nrow(infile)]
  #Get summary information
  summary <- infile[ which(infile$V1=='SUMMARY' | infile$V1=='COMB_TOT_BP' | infile$V1=='COMB_TOT_PROP'), ]
    if (nrow(summary) < 2) stop ("SUMMARY information not present, there is something wrong with the input file")
  summary <- summary[,-1]
  namez <- paste0(outname, "_summary.txt")
  write.table(summary, namez , sep = " ", col.names = F, row.names = F, quote =F)
  summary <- as.data.frame(summary)
  print("Summary Written")
  
  #Get genome sizes in basepairs
  #will get NAs, don't worry
  gsize = as.numeric(as.character(summary[1:2,5]))
  asize =gsize[2]
  gsize = gsize[1]

  #make rounding function
  mround <- function(x,base){ 
    base*round(x/base) 
  } 

  #make a variety of color ramp for  graphs
  colfunc1 <- colorRampPalette(c("darkred", "darkgreen"))
  colfunc2 <- colorRampPalette(c("darkred", "darkblue"))
  colfunc3 <- colorRampPalette(c("orange", "darkblue"))
  colfunc4 <- colorRampPalette(c("darkgreen", "red", "darkblue", "orange", "purple", "yellow", "blue", "red", "green"))
  colfunc5 <- colorRampPalette(c("gold", "darkblue", "red", "orange", "green"))
  colfunc6 <- colorRampPalette(c("green", "darkgreen"))
  colfunc7 <- colorRampPalette(c("red", "darkred"))
  colfunc8 <- colorRampPalette(c("blue", "darkblue"))
  colfunc9 <- colorRampPalette(c("yellow", "orange"))
  colfunc10 <- colorRampPalette(c("white", "black"))
  colfunc11 <- colorRampPalette(c("white", "gray", "black"))

  #Total mean coverage for populations 
  tma <-   infile[ which(infile$V1=='TMA'), ]
  tma<- tma[1:3]
  tms <-  infile[ which(infile$V1=='TMS'), ]
  tma$V4 <- tms$V3
  tma <- tma[,-1]
  colnames(tma) <- c("Population", "Mean_Coverage", "Stdev")
  tma$Mean_Coverage <- as.numeric(as.character(tma$Mean_Coverage))
  tma$Stdev <- as.numeric(as.character(tma$Stdev))
  tma$Population <- as.numeric(as.character(tma$Population))
  tma <- tma[order(tma$Population),]
  
    if (nrow(tma) < 1) stop ("Mean coverage information not present; something is wrong with input file")
  
  if (nrow(tma) > 1) {
    meanr <- mean(tma$Mean_Coverage)
    stdr <- stdev(tma$Mean_Coverage)
    sumr <- data.frame("Mean", meanr,stdr)
    colnames(sumr) <- c("Population", "Mean_Coverage", "Stdev")
    tma <- rbind(tma,sumr)
  }
  
  pnumz <- as.numeric(nrow(tms))
  tma$Population <- factor(tma$Population, levels = tma$Population)
  
  #Get better ticks - get ticks based on scale 
  scalehigh <- colMaxs(tma$Mean_Coverage)+colMaxs(tma$Stdev)
  scalelow <- colMins(tma$Mean_Coverage)-colMaxs(tma$Stdev)
  diff <- scalehigh - scalelow
    
    ##TICK ADJUST##
    if (diff > 0 && diff < .05) {
      bey =.0001
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    
    if (diff > 0.05 && diff < .15) {
      bey =.005
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > .15 && diff < .25) {
      bey =.01
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    
    if (diff > .25 && diff < .5) {
      bey =.025
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    
    if (diff > .5 && diff < 1.05) {
      bey =.05
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 1.05 && diff < 2.5) {
      bey =.15
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 2.5 && diff < 5) {
      bey =.25
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    
    if (diff > 5 && diff < 10) {
      bey =.5
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 10 && diff < 20) {
      bey =1
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 20 && diff < 40) {
      bey =2
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 40 && diff < 80) {
      bey =5
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 80 && diff < 100) {
      bey =5
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 100 && diff < 150) {
      bey =10
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 150 && diff < 200) {
      bey =15
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    
    if (diff > 200 && diff < 400) {
      bey =50
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 400) {
      bey =100
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
     ##TICK ADJUST DONE##

    #plot
    q <- ggplot(tma, aes(x = Population, y = Mean_Coverage, fill = Population)) + 
      geom_bar(stat = "identity", width =.9, colour ="black", size =.25)  +
       geom_errorbar(aes(ymin=Mean_Coverage-Stdev,ymax=Mean_Coverage+Stdev), width=.2, size = .25, position=position_dodge(.9))

    p <- q + labs(title= "Mean Depth of Coverage", x= "Library", y = "Mean Depth of Coverage (X)") + 
     theme(
     #legend.position="none",
       axis.ticks.x=element_blank(),  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(),
     panel.border = element_rect(colour = "black", fill=NA, size=.25)) +
       
      scale_fill_manual(values= c(rep("darkblue",pnumz), "darkgreen")) +
     # scale_fill_manual( values= colfunc1(nrow(tma))) +
     geom_hline(yintercept=0)
  
    r <- p + scale_y_continuous(breaks=seq(0,scalez,bey))

    namez <- paste0(outname, "_mean_coverage.pdf")
    savePlot <- function(r) {
      pdf(namez, width=10, height=6)
      print(r)
    
        dev.off()
    }
    savePlot(r)

    print("Mean Coverage Stats Complete")
  
  #Total mean for all populations - Filtered - AKA analysis coverage
  fma <-   infile[ which(infile$V1=='FMA'), ]
  fma<- fma[1:3]
  fms <-  infile[ which(infile$V1=='FMS'), ]
  fma$V4 <- fms$V3
  fma <- fma[,-1]
  colnames(fma) <- c("Population", "Mean_Coverage", "Stdev")
  fma$Mean_Coverage <- as.numeric(as.character(fma$Mean_Coverage))
  fma$Stdev <- as.numeric(as.character(fma$Stdev))
  fma$Population <- as.numeric(as.character(fma$Population))
  fma <- fma[order(fma$Population),]
  
  if (nrow(fma) < 1) stop ("Filtered coverage information not present; something is wrong with input file")
  
  if (nrow(fma) > 1) {
    meanr <- mean(fma$Mean_Coverage)
   stdr <- stdev(fma$Mean_Coverage)
    sumr <- data.frame("Mean", meanr,stdr)
    colnames(sumr) <- c("Population", "Mean_Coverage", "Stdev")
    fma <- rbind(fma,sumr)
  } 
  
  fma$Population <- factor(fma$Population, levels = fma$Population)

 #plot
  scalehigh <- colMaxs(fma$Mean_Coverage)+colMaxs(fma$Stdev)
  scalelow <- colMins(fma$Mean_Coverage)-colMaxs(fma$Stdev)
  diff <- scalehigh - scalelow
  
  ##TICK ADJUST##
  if (diff > 0 && diff < .05) {
    bey =.0001
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > 0.05 && diff < .15) {
    bey =.005
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > .15 && diff < .25) {
    bey =.01
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > .25 && diff < .5) {
    bey =.025
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > .5 && diff < 1.05) {
    bey =.05
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 1.05 && diff < 2.5) {
    bey =.15
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 2.5 && diff < 5) {
    bey =.25
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > 5 && diff < 10) {
    bey =.5
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 10 && diff < 20) {
    bey =1
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 20 && diff < 40) {
    bey =2
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 40 && diff < 80) {
    bey =5
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 80 && diff < 100) {
    bey =5
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 100 && diff < 150) {
    bey =10
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 150 && diff < 200) {
    bey =15
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  
  if (diff > 200 && diff < 400) {
    bey =50
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 400) {
    bey =100
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
 
   ##TICK ADJUST DONE##
  
    q <- ggplot(fma, aes(x = Population, y = Mean_Coverage, fill = Population)) + 
      geom_bar(stat = "identity", width =.9, colour ="black", size =.25) +  
       geom_errorbar(aes(ymin=Mean_Coverage-Stdev,ymax=Mean_Coverage+Stdev), width=.2, size =.25, position=position_dodge(.9))
 
    p <- q + labs(title= "Mean Depth of Coverage After Filters", x= "Population", y = "Mean Depth of Coverage (X)") + 
      theme(
         #legend.position="none",
        axis.ticks.x=element_blank(),  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.25)) +
      scale_fill_manual(values= c(rep("orange",pnumz), "darkgreen")) + 
      #scale_fill_manual( values= colfunc2(nrow(fma))) +
       geom_hline(yintercept=0)
 
    p <- p + scale_y_continuous(breaks=seq(0,scalez,bey))
    
     namez <- paste0(outname, "_mean_filt_cov.pdf")
     savePlot <- function(p) {
       pdf(namez, width=10, height=6)
       print(p)
       dev.off()
     }
     savePlot(p)
     
 #Write out coverage tables
     namez <- paste0(outname, "_mean_cov.txt")
     tma$type <- "No-Filter"
     fma$type <- "Filtered"   
     coma <- rbind(tma,fma)
       write.table(coma, namez , sep = " ", col.names = T, row.names = F, quote =F)

     print("Filtered coverage stats complete")
       
  ################################################      
 #Proportion of genome covered at current filters
 ################################################
       
  fpc <-   infile[ which(infile$V1=='FPC'), ]
    fpc<- fpc[,1:3]
  
  fbc <-  infile[ which(infile$V1=='FBC'), ]
  fpc$V4 <- fbc$V3
  fpc <- fpc[,-1]
  colnames(fpc) <- c("Population", "Proportion_Covered", "Total_BP")
  fpc$Total_BP=0
  fpc$Proportion_Covered <- as.numeric(as.character(fpc$Proportion_Covered))
  fpc$Total_BP <- as.numeric(as.character(fpc$Total_BP))
  fpc$Population <- as.numeric(as.character(fpc$Population))
  fpc <- fpc[order(fpc$Population),]
  
  if (nrow(fpc) < 1) stop ("Proportion information not present; something is wrong with input file")
  
  if (nrow(fpc) > 1){
    meanr <- mean(fpc$Proportion_Covered)
   stdr <- stdev(fpc$Proportion_Covered)
    sumr <- data.frame("Mean", meanr,stdr)
    colnames(sumr) <- c("Population", "Proportion_Covered", "Total_BP")
   fpc <- rbind(fpc,sumr)
  }
 
  fpc$Population <- factor(fpc$Population, levels = fpc$Population)
  
  #plot
  scalehigh <- colMaxs(fpc$Proportion_Covered) +  colMaxs(fpc$Proportion_Covered* .05) 
  scalelow <- 0
  diff <- scalehigh - scalelow
  
  ##TICK ADJUST##
  if (diff > 0 && diff < .05) {
    bey =.0001
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 0.05 && diff < .15) {
    bey =.005
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > .15 && diff < .25) {
    bey =.01
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > .25 && diff < .5) {
    bey =.025
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > .5 && diff < 1.05) {
    bey =.05
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 1.05 && diff < 2.5) {
    bey =.15
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > 2.5 && diff < 5) {
    bey =.25
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 5 && diff < 10) {
    bey =.5
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 10 && diff < 20) {
    bey =1
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > 20 && diff < 40) {
    bey =2
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 40 && diff < 80) {
    bey =5
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 80 && diff < 100) {
    bey =5
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 100 && diff < 150) {
    bey =10
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 150 && diff < 200) {
    bey =15
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 200 && diff < 400) {
    bey =50
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 400) {
    bey =100
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  ##TICK ADJUST DONE##
  
  q <- ggplot(fpc, aes(x = Population, y = Proportion_Covered, fill = Population)) + 
     geom_bar(stat = "identity", width =.9, colour ="black", size =.25) +
     geom_errorbar(aes(ymin=Proportion_Covered-Total_BP,ymax=Proportion_Covered+Total_BP), size = .25, width=.2, position=position_dodge(.9))
  
   p <- q + labs(title= "Proportion of Genome Covered After Coverage Filters", x= "Library", y = "Proportion of Genome") + 
    theme(
     #legend.position="none",
       axis.ticks.x=element_blank(),  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
       panel.border = element_rect(colour = "black", fill=NA, size=.25)) +
   # scale_fill_manual( values= colfunc3(nrow(fpc))) +
     scale_fill_manual(values= c(rep("goldenrod3",pnumz), "gray")) + 
    geom_hline(yintercept=0) 
    # geom_hline(yintercept=.5, linetype =2) #if a cut-off line is desired
  
     t <- p + scale_y_continuous("Proportion of Genome", sec.axis = sec_axis(~.*gsize/1000000, "bp (million)"), breaks=seq(0,scalez,bey))
  
    namez <- paste0(outname, "_prop_cov.pdf")
    savePlot <- function(t) {
      pdf(namez, width=10, height=6)
      print(t)
      dev.off()
    }
    savePlot(t)
    
    
    #Write out coverage tables
    namez <- paste0(outname, "_mean_prop.txt")
    coma <- fpc
    coma[1:pnumz,3]<-as.numeric(as.character(fbc[[3]]))
    coma[(pnumz + 1),3] <- mean(as.numeric(coma[1:pnumz,3]))
    coma <- coma[complete.cases(coma), ]
    write.table(coma, namez , sep = " ", col.names = T, row.names = F, quote =F)
    print ("Proportion of genome stats complete")
    
    
#########################################################
#Anchored vs scaffold proportion of genome after filters  
#########################################################
    
#Proportion of genome covered at current filters
    SCAFP  <-   infile[ which(infile$V1=='SCAFP'), ]
    SCAFP <- SCAFP [,1:4]
    SCAFL  <-  infile[ which(infile$V1=='SCAFL'), ]
    SCAFP$V5 <- SCAFL$V4
    SCAFP <- SCAFP[,-1]
    colnames(SCAFP) <- c("Population", "Region", "Proportion_Covered", "Total_BP")
    SCAFP$Total_BP=0
    SCAFP$Total_BP <- as.numeric(as.character(SCAFP$Total_BP))
    SCAFP$Population <- as.numeric(as.character(SCAFP$Population))
    SCAFP$Proportion_Covered <- as.numeric(as.character(SCAFP$Proportion_Covered))
    SCAFP <- SCAFP[order(SCAFP$Population),]
    
    if (nrow(SCAFP) > 1){
      meanr <- mean(SCAFP$Proportion_Covered)
      stdr <- stdev(SCAFP$Proportion_Covered)
      sumr <- data.frame("13371", "scaff",meanr,stdr)
      colnames(sumr) <- c("Population", "Region", "Proportion_Covered", "Total_BP")
      SCAFP <- rbind(SCAFP,sumr)
    }
    
    ANCP   <-   infile[ which(infile$V1=='ANCP'), ]
    if (nrow(ANCP) > 1){
      ANCP  <- ANCP  [, colSums(ANCP   != "") != 0]
    }
    

    if (nrow(ANCP ) < 2){
      ANCP <- ANCP[,2:5]
      ANCL   <-  infile[ which(infile$V1=='ANCL'), ]
      ANCP$V5 = ANCL[1,4]
      colnames(ANCP ) <- c("Population", "Region", "Proportion_Covered", "Total_BP")
      SCAFP$Total_BP <- SCAFL[1,4]
      ANCP$Proportion_Covered <- as.numeric(as.character(ANCP$Proportion_Covered))
      ANCP$Total_BP<-as.numeric(as.character(ANCP$Total_BP))
      SCAFP$Proportion_Covered <- as.numeric(as.character( SCAFP$Proportion_Covered))
      SCAFP$Total_BP<-as.numeric(as.character(SCAFP$Total_BP))
    }
    
    if (nrow(ANCP) > 1){
      ANCL   <-  infile[ which(infile$V1=='ANCL'), ]
      ANCP$V5 <- ANCL$V4
      ANCP <- ANCP [,-1]
      colnames(ANCP ) <- c("Population", "Region", "Proportion_Covered", "Total_BP")
      ANCP$Total_BP=0
      ANCP$Total_BP <- as.numeric(as.character(ANCP$Total_BP))
      ANCP$Population <- as.numeric(as.character(ANCP$Population))
      ANCP$Proportion_Covered <- as.numeric(as.character(ANCP$Proportion_Covered))
      ANCP <- ANCP[order(ANCP$Population),]
      meanr <- mean(ANCP$Proportion_Covered)
      stdr <- stdev(ANCP$Proportion_Covered)
      sumr <- data.frame("13371", "anch",meanr,stdr)
      colnames(sumr) <- c("Population", "Region", "Proportion_Covered", "Total_BP")
      ANCP <- rbind(ANCP,sumr)
      ANCP$Proportion_Covered <- as.numeric( ANCP$Proportion_Covered)
      ANCP$Total_BP<-as.numeric(ANCP$Total_BP)
      SCAFP$Proportion_Covered <- as.numeric( SCAFP$Proportion_Covered)
      SCAFP$Total_BP<-as.numeric(SCAFP$Total_BP)
    }

  SCANC <- rbind(ANCP,SCAFP)    

  if (nrow(ANCP) > 1){
    SCANC$Population <- as.numeric(as.character(SCANC$Population))
    SCANC <- SCANC[order(SCANC$Population),]
    SCANC[SCANC==13371] <- "Mean"
    SCANC$Region = gsub("anch", "Anchored", SCANC$Region)
    SCANC$Region = gsub("scaff", "Unanchored", SCANC$Region)
  }

  options(warn=-1)
  SCANC$Population <- factor(SCANC$Population, levels = SCANC$Population)
    
  #plot

    scalehigh <- colMaxs(SCANC$Proportion_Covered) +  colMaxs(SCANC$Proportion_Covered* .05) 
    scalelow <- 0
    diff <- scalehigh - scalelow
    
    ##TICK ADJUST##
    if (diff > 0 && diff < .05) {
      bey =.0001
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 0.05 && diff < .15) {
      bey =.005
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > .15 && diff < .25) {
      bey =.01
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > .25 && diff < .5) {
      bey =.025
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > .5 && diff < 1.05) {
      bey =.05
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 1.05 && diff < 2.5) {
      bey =.15
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 2.5 && diff < 5) {
      bey =.25
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    
    if (diff > 5 && diff < 10) {
      bey =.5
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 10 && diff < 20) {
      bey =1
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 20 && diff < 40) {
      bey =2
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 40 && diff < 80) {
      bey =5
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 80 && diff < 100) {
      bey =5
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 100 && diff < 150) {
      bey =10
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 150 && diff < 200) {
      bey =15
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    
    if (diff > 200 && diff < 400) {
      bey =50
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    
    if (diff > 400) {
      bey =100
      scalez <- scaley <- mround(scalehigh,bey)
      scaley <- mround(scalez, bey)
    }
    ##TICK ADJUST DONE##

    q <- ggplot(SCANC, aes(Population,Proportion_Covered, fill = Population)) + 
      geom_bar(aes(fill = Region), position = "dodge", stat = "identity", width =.9, colour ="black", size =.5) +
      scale_fill_manual(values = c("darkslateblue","aquamarine3")) 

    p <- q + labs(title= "Proportion of Genome Covered After Coverage Filters: Anchored Vs. Unanchored Alignment", x= "Library", y = "Proportion of Genome", fill = "Assembly Type") + 
       theme(
        #legend.position="none",
        axis.ticks.x=element_blank(),  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.1)) +
      geom_hline(yintercept=0) +
      geom_hline(yintercept=0.5, linetype = 2) #if cut off line desired
 
    t <- p + scale_y_continuous("Proportion of Genome", sec.axis = sec_axis(~.*gsize, "basepairs"), breaks=seq(0,scalez,bey))
    
    namez <- paste0(outname, "_prop_anc_vs_scaff.pdf")
    savePlot <- function(t) {
      pdf(namez, width=10, height=6)
      print(t)
      dev.off()
    }
    savePlot(t)
    print ("Anchored vs. Scaffold stats complete")
  
####################################################    
#Proportion of genome covered at different depths of coverage
###############################################
    
    dop <-   infile[ which(infile$V1=='DOP'), ]
    dop <- dop[,1:4]
    dop$V2 <- as.numeric(as.character(dop$V2))
    dop$V3 <- as.numeric(as.character(dop$V3))
    dop <- dop[order(dop$V3,dop$V2),]
    ctab <- as.data.frame (unique(dop$V2))
    colnames(ctab) <- "Coverage"
    
    if (nrow(dop) < 1) stop ("Proportion information not present; something is wrong with input file")

    iter=2
    for (i in (unique(dop$V3))) {
      c1 <- dop[ which(dop$V3==i), ]
      c2 <- c1$V4
      ctab[,iter] <- c2
      ctab[,iter] <- as.numeric(as.character(ctab[,iter]))
      colnames(ctab)[iter] <- i
      iter=iter+1
    }

  if (i > 1){
   ctab$mean <- rowMeans(ctab[, -1])
    ctab$stdev <- rowStdevs(ctab[, -1])
   ctab2 <- ctab[,-(ncol(ctab))]
  }
  if (i < 2){  
   ctab2 <- ctab
  }
    
  df.melted <- melt(ctab2, id = "Coverage")
  df.melted$Coverage <- as.numeric(as.character(df.melted$Coverage))
  colnames(df.melted) <- c("Coverage", "Library", "value")

  scalehigh <- 1
  scalelow <- 0
  diff <- scalehigh - scalelow
  ##TICK ADJUST##
  if (diff > 0 && diff < .05) {
    bey =.0001
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > 0.05 && diff < .15) {
    bey =.005
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > .15 && diff < .25) {
    bey =.01
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > .25 && diff < .5) {
    bey =.025
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > .5 && diff < 1.05) {
    bey =.05
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > 1.05 && diff < 2.5) {
    bey =.15
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 2.5 && diff < 5) {
    bey =.25
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > 5 && diff < 10) {
    bey =.5
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > 10 && diff < 20) {
    bey =1
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > 20 && diff < 40) {
    bey =2
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > 40 && diff < 80) {
    bey =5
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > 80 && diff < 100) {
    bey =5
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > 100 && diff < 150) {
    bey =10
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 150 && diff < 200) {
    bey =15
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  if (diff > 200 && diff < 400) {
    bey =50
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  
  if (diff > 400) {
    bey =100
    scalez <- scaley <- mround(scalehigh,bey)
    scaley <- mround(scalez, bey)
  }
  ##TICK ADJUST DONE##

  p <- ggplot(df.melted, aes(x = Coverage, y = value, colour = Library)) +  
    labs(title= "Proportion of Genome Covered at Different Depths") +
    geom_line(size=.25, linetype =1) + 
    ylab(label="Proportion of genome covered") + 
    xlab("Coverage (X)") + 
    scale_colour_manual(values= c(colfunc4(nrow(tma)-1), "black")) +
    scale_linetype_manual(values=c(rep("dashed",(nrow(tma)-1)), "dotted"))
    
 # Dotted line
  #geom_vline(xintercept=seq(0,150,5), color ="grey",  linetype = "longdash", lwd=.25) +
  
  f<- p + scale_y_continuous("Proportion of Genome", sec.axis = sec_axis(~.*gsize/1000000, "bp (million)"), limits = c(0,1), expand = c(0, 0),
     breaks=seq(0,scalez,bey)) +
     scale_x_continuous(breaks=seq(0,150,5), expand = c(0, 0))
  
  f <- f + theme(
    axis.text.x  = element_text(angle=90, vjust=0.5),
   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=.25)) +
    scale_fill_manual( values= colfunc3(nrow(fpc)))
  
  #xaxis ticks
  #axis.ticks.x=element_blank()
           
  namez <- paste0(outname, "_prop_at_covs.pdf")
  savePlot <- function(f) {
    pdf(namez, width=7, height=5)
    print(f)
    dev.off()
  }
  savePlot(f)

  
#Write table out
  dop <-dop[,-1]
  colnames(dop) <- c("Coverage(X)", "Population", "Proportion")  
  namez <- paste0(outname, "_prop_covs.txt")
  write.table(dop, namez , sep = " ", col.names = T, row.names = F, quote =F)
  print("Proportion versus coverage stats complete")
  
####################  
#Chromosome Analyses
##################

 #Get chromosome sizes from summmary
 chrsize <-   infile[ which(infile$V1=='CHROMOSOME'), ]
 scaffsize <-   infile[ which(infile$V1=='SCAFFOLD'), ]
  if (nrow(chrsize) < 2) stop ("Chromosome stats not present; will skip chromosome analyses")
 s2 <- as.numeric(as.character(scaffsize[1,2]))
 sciff <- data.frame("scaff", s2)       
 colnames(sciff) <- c('chr', 'bp')
 chrsize <- data.frame(do.call('rbind', strsplit(as.character(chrsize$V2),'\t',fixed=TRUE)))
 colnames(chrsize) <- c('chr', 'bp')
 chrsize$bp <- as.numeric(as.character(chrsize$bp))
 chrsize <- rbind(chrsize, sciff)

 #Propotion of alignments between coverage filters 
 chrp <-   infile[ which(infile$V1=='CHRP'), ]
 chrp<- chrp[, 1:4]
 mchr <-   infile[ which(infile$V1=='MCHR'), ]
 mchr<- mchr[, 1:4]
 schr <-   infile[ which(infile$V1=='SCHR'), ]
 mchr$V5 <- schr$V4
 mchr <- mchr[,-1]
 colnames(mchr) <- c("Population", "Chromosome", "Mean_Cov", "stdev")
 
 #..and now scaffolds
 scafz <-   infile[ which(infile$V1=='SCOV'), ]
 scafz<- scafz[, 1:4]
 scafd <- infile[ which(infile$V1=='SDEV'), ]
 scafz$V5<- scafd$V4
 scafz <- scafz[,-1]
 colnames(scafz) <- c("Population", "Chromosome", "Mean_Cov", "stdev")
 
 #combine scaffold info with anchored info
 mchr <- rbind(scafz, mchr)
 
  ###############################
 ### CHROMOSOME PROPORTIONS ####
 ###############################
 
 #Get chromosome sizes from summmary
 chrsize <-   infile[ which(infile$V1=='CHROMOSOME'), ]
 scaffsize <-   infile[ which(infile$V1=='SCAFFOLD'), ]
 s2 <- as.numeric(as.character(scaffsize[1,2]))
 sciff <- data.frame("scaff", s2)       
 colnames(sciff) <- c('chr', 'bp')
 chrsize <- data.frame(do.call('rbind', strsplit(as.character(chrsize$V2),'\t',fixed=TRUE)))
 colnames(chrsize) <- c('chr', 'bp')
 chrsize$bp <- as.numeric(as.character(chrsize$bp))
 chrsize <- rbind(chrsize, sciff)
 
 #Propotion of alignments between coverage filters 
 chrp <-   infile[ which(infile$V1=='CHRP'), ]
 chrp<- chrp[,2:4]
 colnames(chrp) <- c("Population", "Chromosome", "Proportion")
 
 #..and now scaffolds
 scafz <-   infile[ which(infile$V1=='SCAFP'), ]
 scafz<- scafz[, 2:4]
 colnames(scafz) <- c("Population", "Chromosome", "Proportion")
 
 #combine scaffold info with anchored info
 mchr <- rbind(scafz, chrp)
 
###MEANS ONLY
 
 #Get mean values for all populations first
 mchr$Proportion <- as.numeric(as.character(mchr$Proportion))
 mean.mchr <- aggregate(mchr[, 3], list(mchr$Chromosome), mean)
 std.mchr <- aggregate(mchr[, 3], list(mchr$Chromosome), stdev)
 mean.mchr$chr <- chrsize$chr
 mean.mchr$mean <- mean.mchr$x
 mean.mchr$stdev <- std.mchr$x
 mean.mchr$size <- chrsize$bp
 mean.mchr<- mean.mchr[,-(1:2)]
 mean.mchr$size <- as.numeric(as.character(mean.mchr$size))
 mean.mchr$prop <- mean.mchr$size / sum(mean.mchr$size) *10
 
 #plot means
 mean.mchr$w <- cumsum(mean.mchr$prop)
 mean.mchr$wm <- mean.mchr$w - mean.mchr$prop
 mean.mchr$wt <- with(mean.mchr, wm + (w - wm)/2)
 mean.mchr$w2 <- cumsum(mean.mchr$prop)
 mean.mchr$wm2 <- mean.mchr$w - mean.mchr$prop
 mean.mchr$wt <- with(mean.mchr, wm + (w - wm)/2)
 
 scalehigh <- colMaxs(mean.mchr$mean)+ colMaxs(mean.mchr$mean* .5) 
 scalelow <- 0
 diff <- scalehigh - scalelow
 ##TICK ADJUST##
 if (diff > 0 && diff < .05) {
   bey =.0001
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 if (diff > 0.05 && diff < .15) {
   bey =.005
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 if (diff > .15 && diff < .25) {
   bey =.01
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 if (diff > .25 && diff < .5) {
   bey =.025
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 if (diff > .5 && diff < 1.05) {
   bey =.05
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 if (diff > 1.05 && diff < 2.5) {
   bey =.15
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 if (diff > 2.5 && diff < 5) {
   bey =.25
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 if (diff > 5 && diff < 10) {
   bey =.5
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 if (diff > 10 && diff < 20) {
   bey =1
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 if (diff > 20 && diff < 40) {
   bey =2
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 if (diff > 40 && diff < 80) {
   bey =5
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 if (diff > 80 && diff < 100) {
   bey =5
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 if (diff > 100 && diff < 150) {
   bey =10
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 if (diff > 150 && diff < 200) {
   bey =15
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 if (diff > 200 && diff < 400) {
   bey =50
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 
 if (diff > 400) {
   bey =100
   scalez <- scaley <- mround(scalehigh,bey)
   scaley <- mround(scalez, bey)
 }
 
 ##TICK ADJUST DONE##
 
  p <- ggplot(mean.mchr, aes(ymin = 0))
  q <- p + geom_rect(aes(xmin = wm, xmax = w, fill =chr, ymax = mean), colour="black")
  p1 <- q  +  scale_fill_manual(values= c(rep("grey",(nrow(mean.mchr)))))
  p2 <- p1 + geom_text(aes(x = wt, y = mean *  0.5, label = chr, angle = 90, hjust = 0.5, size=.75)) 
  p3 <- p2 + theme_bw() + theme(legend.position = "none",  axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
                               panel.border = element_rect(colour = "black", fill=NA, size=.25)) +
   labs(title = "Proportion of each chromosome covered at across all libraries", y= "Mean Coverage", x="Assembly Chromosome" ) +
   geom_hline(yintercept= mean(mean.mchr$mean), size=1, linetype ="dashed", color ="darkred" ) +  geom_hline(yintercept= 0 ) 
 
 if (i >1 ) {
   f <-p3 + geom_errorbar(data =mean.mchr, aes(ymin = mean -stdev, ymax = mean + stdev, x = wt, width= min(w*.10)),
                          position = position_dodge(width = .9)) 
 }
 if (i < 2) {
   f <- p3  
 }
 f <- f + scale_y_continuous("Proportion of Genome", sec.axis = sec_axis(~.*gsize/1000000, "bp (million)"), breaks=seq(0,scalez,bey)) + scale_x_continuous( expand = c(0, 0))
 
 namez <- paste0(outname, "_chromosome_prop_Mean.pdf")
 savePlot <- function(f) {
   pdf(namez, width=16, height=6)
   print(f)
   dev.off()
 }
 savePlot(f)
 
 #Write table out
 chrinfo <-  mean.mchr[,1:4]
 namez <- paste0(outname, "_chromosome_prop_Mean.txt")
 write.table(chrinfo, namez , sep = " ", col.names = T, row.names = F, quote =F)
   
  print("Chromosome mean stats complete")
 
 ############################
 #Now, for each population
 ############################
 popz <- as.vector(unique(mchr$Population))
 
 if (length(popz) >1) {
 pI=0
 for (i in (popz)) {
   pI = pI + 1
   subz <-   mchr[ which(mchr$Population==i), ]
   subz$order <- as.numeric(gsub("[^0-9]", "", subz$Chromosome)) 
   subz <- subz[order(subz$order),] 
   subz$size <- chrsize$bp
   subz$prop <- subz$size / sum(subz$size) *10
   subz$w <- cumsum(subz$prop)
   subz$wm <- subz$w - mean.mchr$prop
   subz$wt <- with(subz, wm + (w - wm)/2)
   subz$w2 <- cumsum(subz$prop)
   subz$wm2 <- subz$w - subz$prop
   subz$wt <- with(subz, wm + (w - wm)/2)
   colorz <- colfunc4(length(popz))
   c1 <- (colorz[pI])
   colorz2 <- replicate(nrow(subz), c1)   
   
   scalehigh <- colMaxs(subz$Proportion)+ colMaxs(subz$Proportion) *.95 
   scalelow <- 0
   diff <- scalehigh - scalelow
   ##TICK ADJUST##
   if (diff > 0 && diff < .05) {
     bey =.0001
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   
   
   if (diff > 0.05 && diff < .15) {
     bey =.005
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   
   if (diff > .15 && diff < .25) {
     bey =.01
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   
   
   if (diff > .25 && diff < .5) {
     bey =.025
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   
   
   if (diff > .5 && diff < 1.05) {
     bey =.05
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   
   if (diff > 1.05 && diff < 2.5) {
     bey =.15
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   
   if (diff > 2.5 && diff < 5) {
     bey =.25
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   
   
   if (diff > 5 && diff < 10) {
     bey =.5
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   
   if (diff > 10 && diff < 20) {
     bey =1
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   
   if (diff > 20 && diff < 40) {
     bey =2
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   
   if (diff > 40 && diff < 80) {
     bey =5
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   
   if (diff > 80 && diff < 100) {
     bey =5
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   
   if (diff > 100 && diff < 150) {
     bey =10
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   
   if (diff > 150 && diff < 200) {
     bey =15
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
  if (diff > 200 && diff < 400) {
     bey =50
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   
   if (diff > 400) {
     bey =100
     scalez <- scaley <- mround(scalehigh,bey)
     scaley <- mround(scalez, bey)
   }
   ##TICK ADJUST DONE##
   
   Title = paste0("Proportion of Coverage of Population #", i)
   p <- ggplot(subz, aes(ymin = 0))
   q <- p + geom_rect(aes(xmin = wm, xmax = w, fill =Chromosome, ymax = Proportion), colour="black")
   
   p1 <- q  +  scale_fill_manual(values= colorz2)
   
   p2 <- p1 + geom_text(aes(x = wt, y = Proportion *  0.5, label = Chromosome, angle = 90, hjust = 0.5, size=1)) 
   
   p3 <- p2 + theme_bw() + theme(legend.position = "none",  axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
                                 panel.border = element_rect(colour = "black", fill=NA, size=.75)) +
     labs(title = Title, y= "Proportion Covered", x="Assembly Chromosome" ) + geom_hline(yintercept= 0) 
   
   f <- p3  
   
   f <- f + scale_y_continuous("Proportion of Genome", sec.axis = sec_axis(~.*gsize, "basepairs"), breaks=seq(0,scalez,bey)) + scale_x_continuous( expand = c(0, 0))
   

   namez <- paste0(outname, "_chromosome_prop_Pop_", i,".pdf")
   savePlot <- function(f) {
     pdf(namez, width=16, height=6)
     print(f)
     dev.off()
   }
   savePlot(f)
  }
 }
  print("Population Chromosome stats complete")
  print("Finished")
  setwd(file.path(currentD))
}
 
 