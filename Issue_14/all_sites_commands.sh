#!/bin/bash
# SMP 67e892d
# QIIME d4817f1
#

SMP_DIR=/Users/yoshikivazquezbaeza/git_sw/smp_clean/
MAPF=${SMP_DIR}StudentMicrobiomeProject-map.txt

echo "gunzipping ..."
gunzip -c unweighted_unifrac_dm.txt.gz > unweighted_unifrac_dm.txt
gunzip -c weighted_unifrac_dm.txt.gz > weighted_unifrac_dm.txt

echo "creating pc files ..."
principal_coordinates.py -i unweighted_unifrac_dm.txt -o unweighted_unifrac_pc.txt
principal_coordinates.py -i weighted_unifrac_dm.txt -o weighted_unifrac_pc.txt

echo "WeeksSinceStart vectors ..."
make_3d_plots.py -i unweighted_unifrac_pc.txt -m ${MAPF} -o all_sites_unweighted_ordered_by_weekssincestart_vectors --custom_axes='WeeksSinceStart' --colorby='AnyTimeseries,SiteID,PersonalID,Age,University,IBD,Pregnant,BodySite,SmokeCigarettes,ArtificialTanning,DrinkAlcohol,WeekNumber,WeeksSinceStart' --add_vectors='SiteID,WeeksSinceStart'
make_3d_plots.py -i weighted_unifrac_pc.txt -m ${MAPF} -o all_sites_weighted_ordered_by_weekssincestart_vectors --custom_axes='WeeksSinceStart' --colorby='AnyTimeseries,SiteID,PersonalID,Age,University,IBD,Pregnant,BodySite,SmokeCigarettes,ArtificialTanning,DrinkAlcohol,WeekNumber,WeeksSinceStart' --add_vectors='SiteID,WeeksSinceStart'
