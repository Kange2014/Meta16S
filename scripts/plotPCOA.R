#! Plot Mothur's PCoA Data with an optional design file for groups in 2d
#!
#! input: axes: the axes file from mothur's pcoa function.
#!		  loading: the loading file from mothur's pcoa function.
#!		  design: the design file for the groups you which to color the points by. defaults to no design file.  
#!		  ax1, ax2: the axes you wish to plot. Defaults to the first 2.
#!		  Lab1, Lab2: titles for the 1st/x and 2nd/y axes.
#!		  psize: Size of the points plotted.

# need to modify xlab and ylab to contain percent value

plotPCOA <- function(axes,loading,design=FALSE, ax1=1, ax2=2, Lab1="PC1", Lab2= "PC2", main="PCoA plot in 2D", psize=I(4.5)){
  require(ggplot2) 
  axi <- read.table(axes,header = T)
  loadings <- read.table(loading,header=T)
  
  Lab1 = paste(Lab1," (",loadings[1,2],"%)",sep="")
  Lab2 = paste(Lab2," (",loadings[2,2],"%)",sep="")

  sample = axi[,1]
  
  if (design==FALSE){
    qplot(data= axi,x= axis1, y=axis2, xlab=Lab1, ylab=Lab2, color = sample, main = main, size = I(psize))
  }
  else {
    design         <- suppressWarnings(read.table(design, header=F))
    tDesign        <- t(design)
    taxi           <- cbind(tDesign[,2],axi)
    names(taxi)[1] <- "partition" 
    groups         <- taxi[1]
    qplot(data=taxi, x=axis1, y=axis2, color=partition, size = I(psize))
  }
  
  ggsave(file="PCoA_plot.png")
}

args<-commandArgs(T)
plotPCOA(axes=args[1],loading=args[2])
