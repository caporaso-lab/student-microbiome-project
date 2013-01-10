To assess variability in alpha diversity over time, I am using the coefficient of variation (CV).
There are directories for each body site which contain both figures and statistical results.
The influence of gender and university on variability was assessed within each body site.
For gender, I ran ttests to test significance.
For University, ANOVA and Tukeys tests were used.
File names that contain # indicate some level of statistical significance (p>0.05).
Mean values for each metric were also calculated and tested for University and Gender differences.
Black bars in each plot are the mean.
The R codes used are below.

##File input##
map=read.table("FILE_NAME.txt", header=TRUE, sep="\t", check.names=FALSE)

##Creates bar for mean in plot##
stat_sum_single <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.y=fun, colour="black", geom=geom, width = 0.4, ymin=-1000,ymax=-1000, ...)
}

library("ggplot2")

##Generates dotplots##
d=ggplot(map, aes(x= Gender, y=CV_observed_species, color=Gender, title = "Forehead - OTUs Observed"))
d=d + geom_point(position=position_jitter(w=0.2), size=3)
d=d + theme_bw()
d=d + ylab("Coefficient of Variation")
d=d + theme(plot.title = element_text(size = rel(2.2)), axis.title.x = element_text(size=rel(1.5)), axis.title.y = element_text(size=rel(1.5)), axis.text = element_text(size=rel(1.2)))
d + stat_sum_single(mean)

##ttest##
ttest=t.test(CV_observed_species ~ Gender, data=map, na.remove=TRUE)
capture.output(ttest, file="Forehead_CV_richness_gender_ttest.txt")

##ANOVA and Tukey's##
AOVModel=aov(CV_observed_species ~ University, data=map)
Anova=summary(AOVModel)
capture.output(Anova, file="ANOVA_Forehead_CV_richness_University.txt")

comparisons=(TukeyHSD(AOVModel))
capture.output(comparisons, file="ANOVA_Forehead_CV_richness_University_Tukeys.txt")