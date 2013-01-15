Boxplots of UniFrac distances within each individual colored by gender (pink = female, blue = male).
Individuals are sorted by median - highest to lowest
Red line is median.
Black line is mean.
R code is below.


# script to compare w/in individual distances among individuals via boxplots
# limit the samples by those existing in the mapping file. the distance
# matrix can have extra and order doesn't matter

compare_within_dists_among_indivs.f = function(dmat_fp,map_fp,indivs_col_name,sort_method,colorBy){
  dmat = read.table(dmat_fp,header=TRUE,sep="\t",row.names=1,check.names=FALSE)
  map = read.table(map_fp,header=TRUE,sep="\t",row.names=1,check.names=FALSE,comment.char="")
  # get subset of dmat in mapping file
  distToKeep=match(row.names(map),row.names(dmat))
  dmat.sub=dmat[distToKeep,distToKeep]
  dmat.sub=as.dist(dmat.sub)
  # convert dmat to 3 column format
  dmat.clms=data.frame(t(combn(labels(dmat.sub),2)),as.numeric(dmat.sub))
  names(dmat.clms)=c("c1","c2","dists")
  # add column with "within"/"between" cat
  c1cat=map[match(dmat.clms[,1],row.names(map)),indivs_col_name]
  c2cat=map[match(dmat.clms[,2],row.names(map)),indivs_col_name]
  cat=seq(from=1,to=length(c1cat))
  cat[c1cat==c2cat]="within"
  cat[c1cat!=c2cat]="between"
  distsWcats=cbind(dmat.clms,cat)
  # subselect only within distances
  distsWcats.wIn=distsWcats[distsWcats[,4]=="within",]
  # add column with individual ID
  indivIDs=map[match(distsWcats.wIn[,1],row.names(map)),indivs_col_name]
  distsWcats.wIn=cbind(distsWcats.wIn,indivIDs)
  # sort plots by sort method
  bySorting=with(distsWcats.wIn,reorder(indivIDs,dists,sort_method))
  # color boxes by metadata
  posColors=c('lightpink2','lightblue2','gray','green','orange')
  factorLevels=map[match(levels(bySorting),map[,'PersonalID']),colorBy]
  colors=as.character(factorLevels)
  for(i in 1:length(levels(factorLevels))){
    colors[colors==levels(factorLevels)[i]]=posColors[i]
    print(paste(levels(factorLevels)[i],": ",posColors[i]))
  }
  # make boxplots
  pdf('within_dists_per_individual.pdf',onefile=TRUE,width=8.5,11)
  boxplot(dists~bySorting,data=distsWcats.wIn,las=1,xlab="Unweighted UniFrac Distance",col=colors,border="black",cex=.75,boxwex=.75,horizontal=TRUE,cex.axis=0.6,main="Forehead")
  abline(v=median(distsWcats.wIn[,'dists']),lty='solid',lwd=2.5,col='red')
  abline(v=mean(distsWcats.wIn[,'dists']),lty='solid',lwd=2.5,col='black')
  dev.off()
}

# Example usage:
dmat_fp = '/Volumes/fiererlab/users/leffj/student_microbiome/CLEAN_u_unifrac_Gut.txt'
map_fp = '/Volumes/fiererlab/users/leffj/student_microbiome/Gut_mega_map.txt'
indivs_col_name = 'SiteID'
sort_method = 'median'
colorBy = 'Gender' # up to 3 levels

compare_within_dists_among_indivs(dmat_fp,map_fp,indivs_col_name,sort_method,colorBy)
