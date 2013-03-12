#!/bin/bash
# SMP 67e892d
# QIIME d4817f1
#

SMP_DIR=/Users/yoshikivazquezbaeza/git_sw/smp/
MAPF=${SMP_DIR}StudentMicrobiomeProject-map.txt
BETADIV="/Users/yoshikivazquezbaeza/git_sw/smp/beta-diversity/"

# echo "gunzipping ..."
gunzip -c ${BETADIV}unweighted_unifrac_dm.txt.gz > ./unweighted_unifrac_dm.txt
gunzip -c ${BETADIV}weighted_unifrac_dm.txt.gz > ./weighted_unifrac_dm.txt

# run the following code in python
# from qiime.filter import sample_ids_from_metadata_description as s
# mapping_file = open('/Users/yoshikivazquezbaeza/git_sw/smp/StudentMicrobiomeProject-map.txt', 'U').readlines()
# sids_tongue = s(mapping_file,'TongueTimeseries:Yes')
# sids_forhead = s(mapping_file,'ForeheadTimeseries:Yes')
# sids_palm = s(mapping_file,'PalmTimeseries:Yes')
# sids_gut = s(mapping_file,'GutTimeseries:Yes')

# fd = open('sample_ids_to_keep.txt', 'w')
# # what we care about is the union of all these sets of sample identifiers
# for sid in list(set(sids_tongue) | set(sids_forhead) | set(sids_palm) | set(sids_gut)):
#     fd.write('%s\n' % sid)

# fd.close()

echo "filtering the distance matrices ..."
filter_distance_matrix.py -i unweighted_unifrac_dm.txt -o unweighted_unifrac_dm.tsonly.txt --sample_id_fp=sample_ids_to_keep.txt -m ${MAPF}
filter_distance_matrix.py -i weighted_unifrac_dm.txt -o weighted_unifrac_dm.tsonly.txt --sample_id_fp=sample_ids_to_keep.txt -m ${MAPF}

echo "creating pc files ..."
principal_coordinates.py -i unweighted_unifrac_dm.tsonly.txt -o unweighted_unifrac_pc.txt
principal_coordinates.py -i weighted_unifrac_dm.tsonly.txt -o weighted_unifrac_pc.txt

echo "WeeksSinceStart vectors ..."
make_3d_plots.py -i unweighted_unifrac_pc.txt -m ${MAPF} -o all_sites_ts_only_unweighted_ordered_by_weekssincestart_vectors --custom_axes='WeeksSinceStart' --add_vectors='PersonalIDBodySite,WeeksSinceStart'
make_3d_plots.py -i weighted_unifrac_pc.txt -m ${MAPF} -o all_sites_ts_only_weighted_ordered_by_weekssincestart_vectors --custom_axes='WeeksSinceStart' --add_vectors='PersonalIDBodySite,WeeksSinceStart'
