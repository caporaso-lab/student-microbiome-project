Included in this directory are results correlating beta diversity across body habitats.
Each point represents the mean value within an individual for each respective body habitat.
In the plots, points are colored by gender.
Linear model results are presented with F-test values giving us an indication of how well the model explains the variance in the data.
The R-code for the plots and model is below.

##Plot beta vs beta
c=ggplot(map, aes(x=MEAN_Forehead_Unweighted, y=MEAN_Palm_Unweighted, color=Gender))
c=c + geom_point(size=4)
c=c + theme_bw(base_size=20)
line=coef(lm(MEAN_Palm_Unweighted ~ MEAN_Forehead_Unweighted, data=map))
c=c + geom_abline(intercept=(line[1]), slope=(line[2]))
c=c + labs(title = "Forehead vs. Palm", x= "Average Unweighted UniFrac Distance - Palm", y= "Average Unweighted UniFrac Distance - Tongue")
c


##To get results from model, y=mx+b, m=slope, b=y-intercept
model=lm(MAD_Palm_Unweighted ~ MAD_Forehead_Unweighted, data=map)
x=summary(model)
capture.output(x, file="Forehead_vs_Palm_UUniFrac_model_results.txt")