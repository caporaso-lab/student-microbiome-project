I was interested to see how the community changes through time relative to the previous
time point. This could allow us to identify disturbance events.
The code is below.

#Extract unifrac distances for each individual and plot relative to the previous time point.
ExtractSerialCommunityDistance.f=function(map_fp, dist_fp){
	
	dist=read.table(dist_fp, header=TRUE, row.names=1, sep="\t",check.names=FALSE)
	map=read.table(map_fp, header=TRUE, check.names=FALSE, sep="\t")

	pdf("DistanceFromPreviousTime.pdf", onefile=TRUE)
	u=unique(map[,"PersonalID"])
	
	out=NULL
	
	for(i in 1:length(u)){
	  # pull out and order map file entries specific to individual
		tmp=map[map[,"PersonalID"]==u[i],]
    s=sort(tmp[,'WeeksSinceStart'],index=TRUE)
    tmp.ordered=tmp[s$ix,]
    # get matching distance matrix values
    dist.toKeep = match(tmp.ordered[,"SampleID"],row.names(dist))
		tmp.d=dist[dist.toKeep,dist.toKeep]
    # get corresponding distances from previous time points
    distVals=NULL
    for(k in 1:nrow(tmp.ordered)){
#       if(k==1) {
#         distVals=0
#         next
#       }
      distVal=tmp.d[k,k-1]
      distVals=c(distVals,distVal)
    }
    
    indiv.out=cbind(rep(as.character(u[i]),times=length(distVals)),tmp.ordered[2:nrow(tmp.ordered),'WeeksSinceStart'],distVals)
    colnames(indiv.out)=c("Individual","WeeksSinceStart","distFromPrevious")
		
		out=rbind(out,indiv.out)
		
		plot(x=indiv.out[,2],y=indiv.out[,3], ylim=c(.2,1), xlim=c(0,10), main=u[i], type="b", col="coral", pch=19, xlab="WeeksSinceStart", ylab="Unweighted UniFrac Distance From Previous - Forehead")
	}
	
	write.table(out,"DistFromPrevious.txt", sep="\t", quote=FALSE,row.names=FALSE)
	dev.off()
	}