#! Plot Rarefaction.single data which is the output from mothur's rarefraction.single function.
#!
#! input: file: the .group.rarefaction file which you wish to plot. Only works for the output of rarefaction.single from mothur.
#!		  groups: the specific groups within x which you wish to plot, in list form c("Examplegroup1", "Examplegroup2');
#!				  either the names ie c("m24626d36", "Day3") or the number the group falls in the list of all groups,
#!				  ie c(1,4,7) (the 1st 4th and 7th group).
#!		  error: if error = TRUE then plot the Lower and Upperbounds for each group. 
#!		  delNum: the number of excess characters before each name generally X0.03, 
#!				  or something similiar. Set delNum to the number you which to delete.

rarePlot <- function(file,png,groups,ylab= "Number of Different OTUs",xlab= "Number of Tags Sampled",pch= NA, xlim = NULL, ylim = NULL, error=FALSE) #Same Defaults as matplot
{
  library(reshape2)
  require(ggplot2)#install.packages("ggplot2")
  
  png(file = png,res=600,width=4800, height=4800)
  
  Data <- read.table(file, header=T)
  Length <- length(Data[1,])
  if (error == FALSE)
  {
    Steps <- seq(from=2,to=Length-2,by=3)
    TrueDat <- Data[,c(1,Steps)]
    
	# substr(names(TrueDat), 1, 6) <- ""
    # substr(names(TrueDat), 1, 6) <- ""
    # substr(names(TrueDat), 1, 6) <- ""
    # substr(names(TrueDat), 1, 6) <- ""
    # substr(names(TrueDat), 1, 6) <- ""
    # substr(names(TrueDat), 1, 6) <- ""
	# 
	# it seems that substr doesn't work as expected in the plugin
	
	names(TrueDat) <- gsub("[X\\.0-9]{6}","", names(TrueDat))
    names(TrueDat)[1] <- "numSampled"
    
    ## Plot data with ggplot2

    # transform data to long form
    longRareData <- melt(TrueDat, id.vars = "numSampled")

    # plotting
    options(scipen = 100000)
    p <- ggplot(data = longRareData,
       aes(x = numSampled, y = value, color = variable)) +
       geom_line() + 
       ggtitle("Rarefaction Curves of All Samples\n") +
       labs(x = xlab, y = ylab) +
       guides(color = guide_legend("Legend")) +
       theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    print(p)
  }
  else 
  {
    PureSteps <- seq(from=2,to=Length-2,by=3)
    LowSteps <- seq(from=3,to=Length-1,by=3)
    HighSteps <- seq(from=4,to=Length,by=3)
    PureDat <- Data[,c(1,PureSteps)]
    LowDat <- Data[,c(1,LowSteps)]
    HighDat <- Data[,c(1,HighSteps)]
    PureNames <- names(PureDat)
    # substr(names(PureDat), 1, 6) <- ""
    # substr(names(PureDat), 1, 6) <- ""
    # substr(names(PureDat), 1, 6) <- ""
    # substr(names(PureDat), 1, 6) <- ""
    # substr(names(PureDat), 1, 6) <- ""
    # substr(names(PureDat), 1, 6) <- ""
	names(PureDat) <- gsub("[X\\.0-9]{6}","", names(PureDat))
    names(PureDat)[1] <- "Number"
    PureColumns <- PureDat[-1]
    matplot(Number,PureColumns[,groups],ylab= ylab,xlab = xlab, pch=pch,col=1:15, xlim = xlim, ylim = ylim)
    legend("topleft",legend=names(PureColumns[,groups]), col =1:15,pch=19)
    matlines(Number,PureColumns[,groups], col=1:15)
    LowColumns<- LowDat[-1]
    matlines(Number,LowColumns[,groups], col=1:15, lwd=.2)
    HighColumns<- HighDat[-1]
    matlines(Number,HighColumns[,groups], col=1:15, lwd=.2)
  }
  
  dev.off()
}

args<-commandArgs(T)
rarePlot(file=args[1],png=args[2],ylab=args[3])
