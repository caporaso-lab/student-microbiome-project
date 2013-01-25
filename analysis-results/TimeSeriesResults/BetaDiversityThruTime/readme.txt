Included here are plots of beta through time for each individual and mantel tests results
correlating distance matrices with time. The plots have "Weeks Between Samples" on the
x-axis to test for a distance decay relationship. Both weighted and unweighted results are
presented. Red lines in plots are mean distances for that individual.

The code used to generate the plots is here.

##Plotting Community Distances versus time
plotDist.f=function(dist_fp, time_fp, map_fp){
	dist=read.table(dist_fp, header=TRUE, row.names=1, sep="\t",check.names=FALSE)
	time=read.table(time_fp, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
	map=read.table(map_fp, header=TRUE, check.names=FALSE, sep="\t")

  library(vegan)
	
	u=unique(map[,"PersonalID"])
	
	
	pdf("UUnifrac_Time_plots.pdf", onefile=TRUE)
	for(i in 1:length(u)){
		g=grep(u[i],map[,"PersonalID"])
		dist.toKeep = match(map[,"SampleID"][map[,"PersonalID"]==u[i]],row.names(dist))
		dist.t=dist[dist.toKeep,dist.toKeep]
    
		time.toKeep = match(map[,"SampleID"][map[,"PersonalID"]==u[i]],row.names(time))
    	time.t=time[time.toKeep,time.toKeep]
		
		dist.d=as.dist(dist.t)
		time.d=as.dist(time.t)
		d.avg=mean(dist.d)

		plot(x=time.d, y=dist.d, main=u[i], ylim=c(0,1), col="blue", xlab ="Weeks Between Samples", ylab="Unweighted UniFrac Distance", xlim=c(0, 10))
		abline(h=d.avg, col="red", lty="dashed", lwd=1.5)
		
	}
	dev.off()
	
}


The R code for the mantel tests is here. 

#To test for a correlation of similarity between samples with time between samples: Mantel test between distance matrix and time matrix.  Function loops for each site within a dataset.
DistTimeMantel.f=function(time_fp, dist_fp, map_fp){
	time=read.table(time_fp, header=TRUE, row.names=1, sep="\t",check.names=FALSE)
	dist=read.table(dist_fp, header=TRUE, row.names=1, sep="\t",check.names=FALSE)
	map=read.table(map_fp, header=TRUE, check.names=FALSE, sep="\t")
  
	u=unique(map[,"PersonalID"])
	library(vegan)
	out=NULL
	
	for(i in 1:length(u)){
		g=grep(u[i],map[,"PersonalID"])
		time.toKeep = match(map[,"SampleID"][map[,"PersonalID"]==u[i]],row.names(time))
    	t=time[time.toKeep,time.toKeep]
		dist.toKeep = match(map[,"SampleID"][map[,"PersonalID"]==u[i]],row.names(dist))
		d=dist[dist.toKeep,dist.toKeep]
    
		t.d=as.dist(t)
		d.d=as.dist(d)
		print(u[i])
		m=mantel(t.d,d.d, method="pearson", permutations=1000)
		print(m)
		
		m2=mantel(t.d,d.d, method="spearman", permutations=1000)
		print(m2)
		
		m3=mantel(t.d,d.d, method="kendall", permutations=1000)
		print(m3)
		
		
		#o=c(m$statistic, m$signif)
		o=c(m$statistic, m$signif, m2$statistic, m2$signif, m3$statistic, m3$signif)
		out=rbind(out,o)
		}
		row.names(out)=u
		#colnames(out)=c("SpearmanR", "SpearmanP")
		colnames(out)=c("PearsonR", "PearsonP", "SpearmanR", "SpearmanP", "KendallR", "KendallSig")
		write.table(out, "DistTimeMantel_Unweighted.txt", quote=FALSE, sep="\t")
}

The time distance matrix was generated from the WeeksSinceStart column in the mapping file
using the script below.

makeTimeDist.f=function(map_fp){
	map=read.table(map_fp,header=TRUE, check.names=FALSE, sep="\t")
	temp=as.matrix(map[,"WeeksSinceStart"])
	names(temp)=map[,"SampleID"]
	temp.d=dist(temp, method="manhattan", diag=FALSE)
	head(temp.d)
	temp.out=as.matrix(temp.d)
	colnames(temp.out)=map[,"SampleID"]
	row.names(temp.out)=map[,"SampleID"]
	write.table(temp.out, "CLEAN_TimeDist.txt", quote=FALSE, sep="\t")
}