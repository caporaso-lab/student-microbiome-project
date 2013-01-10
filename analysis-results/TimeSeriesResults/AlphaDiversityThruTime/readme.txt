Figures generated for each individual using three alpha diversity metrics (Shannon, PD, richness).
Each value was determined from 10 replicate drawings of 10,000 sequences per sample and imported into mapping file.
Red line in each plot is mean for that individual across time points.
Only "Time Series" individuals were included in the analysis.
The R code used for these plots is below.

#Plotting alpha diversity in time
plotAlphaDiversityThruTime.f=function(mega_map_fp){
  mega=read.table(mega_map_fp, header=TRUE,sep="\t", check.names=FALSE)
  pdf("AlphaDivThruTime.pdf", onefile=TRUE, width=8, height=4)
  u=unique(mega[,"PersonalID"])
  for(i in 1: length(u)){
    tmp=mega[mega[,"PersonalID"]==u[i],]
    par(mfrow=c(1,3))
    plot(x=tmp[,"WeeksSinceStart"],y=tmp[,"ShannonEvenness"], main=u[i], ylab="Shannon Index", col="blue", cex=1.5, xlab="Time(weeks)",ylim=c(0,max(mega[,"ShannonEvenness"])))
    abline(h=mean(tmp[,"ShannonEvenness"]),lty="dashed", lwd=1.5, col="red")
    
    plot(x=tmp[,"WeeksSinceStart"],y=tmp[,"PD"], main=u[i], ylab="Faith's Phylogenetic Diversity", col="blue", cex=1.5, xlab="Time(weeks)",ylim=c(0, max(mega[,"PD"])))
    abline(h=mean(tmp[,"PD"]),lty="dashed", lwd=1.5, col="red")
    
    plot(x=tmp[,"WeeksSinceStart"],y=tmp[,"ObservedOTUs"], main=u[i], ylab="No. OTUs (richness)", col="blue", cex=1.5, xlab="Time(weeks)",ylim=c(0,max(mega[,"ObservedOTUs"])))
    abline(h=mean(tmp[,"ObservedOTUs"]),lty="dashed", lwd=1.5, col="red")
    
  }	
  dev.off()
} 