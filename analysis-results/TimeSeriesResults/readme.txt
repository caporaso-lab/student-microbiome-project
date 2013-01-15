For most of these analyses, I have determine a variety of summary statistics (i.e., mean, MAD, etc.) for both alpha and beta diversity for each individual.
The script used was originally written by Ashley Shade and later modified by Jon Leff and myself.
Depending on the script, the basic input files are a mapping file and a distance matrix. 
Some also require a time distance matrix constructed from the WeeksSinceStart column in the mapping file.

The script used to summarize alpha diversity metrics is below.

#Summarize metrics of diversity
#mega_map_fp="CLEAN_MEGA_map.txt"
divSum.f=function(mega_map_fp){
	mega=read.table(mega_map_fp, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
	u=unique(mega[,"PersonalID"])
	div.out=NULL
	for(i in 1:length(u)){
		tmp=mega[mega[,"PersonalID"]==u[i],]

		shannon=mean(tmp[,"ShannonEvenness"])
		s.sd=sd(tmp[,"ShannonEvenness"])
		s.v=var(tmp[,"ShannonEvenness"])
    s.med=median(tmp[,"ShannonEvenness"])
    s.CV=s.sd/shannon
		rich=mean(tmp[,"ObservedOTUs"])
		rich.sd=sd(tmp[,"ObservedOTUs"])
		rich.v=var(tmp[,"ObservedOTUs"])
    rich.med=median(tmp[,"ObservedOTUs"])
    rich.CV=rich.sd/rich
		pd=mean(tmp[,"PD"])
		pd.sd=sd(tmp[,"PD"])
		pd.v=var(tmp[,"PD"])
    pd.med=median(tmp[,"PD"])
    pd.CV=pd.sd/pd
			out=c(paste(u[i]),shannon,s.sd,s.v,s.med,s.CV,rich,rich.sd,rich.v,rich.med,rich.CV,pd,pd.sd,pd.v,pd.med,pd.CV)
		div.out=rbind(div.out,out)
		
		#print(head(div.out))
		}
	colnames(div.out)=c("PersonalID","Mean_Shannon", "SD_Shannon", "Var_Shannon", "Median_Shannon", "CV_Shannon", "Mean_observed_species", "SD_observed_species", "Var_observed_species", "Median_observed_species", "CV_observed_species", "Mean_PD_whole_tree", "SD_PD_whole_tree", "Var_PD_whole_tree", "Median_PD_whole_tree", "CV_PD_whole_tree")
	write.table(div.out, "DiversitySummary.txt", quote=FALSE, sep="\t", row.names=FALSE)
}



The script used to summarize within individual beat diversity metrics is below.

#calculate an overall variance, standard deviation, mean, and median absolute deviation (MAD) of the unifrac distance for each biome
varDistance.f=function(dist_fp, map_fp){
	map=read.table(map_fp, header=TRUE, sep="\t", check.names=FALSE)
	##map=read.table("Forehead_Mega_Map_Final.txt", header=TRUE, sep="\t", check.names=FALSE)
	dist=read.table(dist_fp, header=TRUE, row.names=1, sep="\t",check.names=FALSE)
	##dist=read.table("CLEAN_W_unifrac_Forehead.txt", header=TRUE, row.names=1, sep="\t",check.names=FALSE)
  
  library(plotrix)
  
	u=unique(map[,"PersonalID"])
	
	v.out=NULL
	for(i in 1:length(u)){
	##dist.tmp=dist[map[,"PersonalID"]==u[i],map[,"PersonalID"]==u[i]]
    ##dist.toKeep = match(map[,"SampleID"][map[,"SiteID"]==u[i]],row.names(dist))
	dist.toKeep = match(map[,"SampleID"][map[,"PersonalID"]==u[i]],row.names(dist))
	dist.tmp=dist[dist.toKeep,dist.toKeep]
	v=var(as.dist(dist.tmp))
	stdev=sd(as.dist(dist.tmp))
	mad=mad(as.dist(dist.tmp))
  mean=mean(as.dist(dist.tmp))
	
	v2=c(v,stdev,mad,mean)
	
	v.out=rbind(v.out,v2)
	
		}
	row.names(v.out)=u
	colnames(v.out)=c("Var", "SD", "MAD", "MEAN")	
	write.table(v.out, "VarianceTable.txt", sep="\t", quote=FALSE)
	
} 