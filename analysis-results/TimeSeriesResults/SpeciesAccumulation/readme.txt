Species accumulation curves are presented for each body habitat for each individual. In 
each  plot, two methods were used. The points represent the "collector" method and is
derived from the samples in the order they were collected. The line is the 
"random" method and randomly selects samples within an individual with no regard to the 
order they were collected. The line in each plot is the average of 100 permutations. Both
plots and figures are included in the directory. Codes are below.


#SpeciesAccumulationCurvesforEventDetection
sac.f=function(otu_fp, map_fp){
	otu=read.table(otu_fp, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
	map=read.table(map_fp, header=TRUE, sep="\t", check.names=FALSE)
	  
  
  # there should be no extra sample ids in the mapping file that aren't in the otu table!!
  # this means you can select which samples to include by deleting them from the mapping file
  
	library(vegan)
		
	#remove RDP ID
	l=ncol(otu)
	p=otu[,l]
	otu=otu[,-l]
	
	u=unique(map[,"PersonalID"])
	
	#pdf("SAC_log.pdf", onefile=TRUE)

	out=NULL
	for(i in 1:length(u)){
		print(u[i])
    #matches PersonalIDs to SampleIDs (column 1 in mapping file) with names in otu table (row 1)
		tmp=otu[,match(map[map[,"PersonalID"]==u[i],1],names(otu))]
    #converts to presence absence otu table
		tmp.pa=1*(tmp>0)
		map.tmp=map[map[,"PersonalID"]==u[i],]
	
    #sort consecutively		
		s=sort(map.tmp[,"WeeksSinceStart"],index=TRUE)
		tmp=tmp[,s$ix]
			
		##s.c=specaccum(t(tmp), method="exact")
		s.c=specaccum(t(tmp), method="collector")
		s.r=specaccum(t(tmp), method="random")
		
		#c is collector's, r is random
		c=s.c$richness
		r=s.r$richness
		t=s$x
		
		si=NULL
		for(y in 1:length(c)){
			si[y]=print(as.character(u[i]))
		}
	
		o=cbind(si,c,r,t)
		colnames(o)=c("PersonalID", "c", "r", "t")
		out=rbind(out,o)
		}
		#return(out)
		write.table(out, "SAC.txt", quote=FALSE, sep="\t")
}
	
	
	##Plotting SAC	
sac_fp="SAC.txt"
function(sac_fp){
  #plotting all together
  sac=read.table(sac_fp, header=TRUE, sep="\t", check.names=FALSE)
  #sac=read.table("SAC.txt", header=TRUE, sep="\t", check.names=FALSE)
  par(mar=(c(5,5,1,2)+0.1))
  u=unique(sac[,"PersonalID"])
  
  mx.x=max(sac[,"t"])
  mn.x=min(sac[,"t"])
  mx.y=max(sac[,"r"])
  mn.y=min(sac[,"r"])
  
  colors=rainbow(length(u))
  pdf("SAC.pdf", onefile=TRUE)
  
  #allplots
  for(i in 1:length(u)){
    tmp=sac[sac[,"PersonalID"]==u[i],]
    if( i == 1){
      plot(x=tmp$t,y=tmp$c,type="p", xlab="Time points (WeeksSinceStart)", ylab="Cumulative richness (No. OTUs)", col=colors[i], xlim=c(0,mx.x), ylim=c(0,mx.y), main="Species Accumulation - Gut")
      lines(x=tmp$t,y=tmp$r,lty="dashed", col="black")
    }
    else{
      points(x=tmp$t,y=tmp$c, type="p", cex=0.9, col=colors[i])
      lines(x=tmp$t,y=tmp$r, lty="dashed", col="black")
    }			
  }
  legend("bottomright", cex=0.5, ncol=6, legend=u,col=colors,pch=19)
  
  #individual plots
  for(i in 1:length(u)){
    tmp=sac[sac[,"PersonalID"]==u[i],]
    plot(x=tmp$t,y=tmp$c,type="p", pch=19, cex=1.2, xlab="Time points (WeeksSinceStart)", ylim=c(0,1500), xlim=c(0,mx.x), ylab="Cumulative richness (No. OTUs)", col=colors[i], main=u[i])
    lines(x=tmp$t,y=tmp$r,lty="dashed")
  }
  
  dev.off()
}