#!/bin/bash
# SMP 1b9ba2646b39a8afdabf1951bd530d3247f5bd6d
# QIIME d400522cd5771402e5d7fd6a7756bef24968c3c6
#

# the vectors traces depend an extra column I added called SiteId, which cats
# the PersonalID and the BodySite column
SMP_DIR=/Users/yoshikivazquezbaeza/git_sw/smp/
MAPF=${SMP_DIR}StudentMicrobiomeProject-map.txt

echo "gunzipping ..."
gunzip -c even10000_unweighted_unifrac_dm.txt.gz > even10000_unweighted_unifrac_dm.txt
gunzip -c even10000_weighted_unifrac_dm.txt.gz > even10000_weighted_unifrac_dm.txt

echo "creating pc files ..."
principal_coordinates.py -i even10000_unweighted_unifrac_dm.txt -o even10000_unweighted_unifrac_pc.txt
principal_coordinates.py -i even10000_weighted_unifrac_dm.txt -o even10000_weighted_unifrac_pc.txt

echo "WeeksSinceStart vectors ..."
make_3d_plots.py -i even10000_unweighted_unifrac_pc.txt -m ${MAPF} -o all_sites_unweighted_ordered_by_weekssincestart_vectors --custom_axes='WeeksSinceStart' --colorby='SiteID,PersonalID,Age,IBD,Pregnant,BodySite,SmokeCigarettes,ArtificialTanning,DrinkAlcohol,WeekNumber,WeeksSinceStart' --add_vectors='SiteID,WeeksSinceStart'
make_3d_plots.py -i even10000_weighted_unifrac_pc.txt -m ${MAPF} -o all_sites_weighted_ordered_by_weekssincestart_vectors --custom_axes='WeeksSinceStart' --colorby='SiteID,PersonalID,Age,IBD,Pregnant,BodySite,SmokeCigarettes,ArtificialTanning,DrinkAlcohol,WeekNumber,WeeksSinceStart' --add_vectors='SiteID,WeeksSinceStart'
