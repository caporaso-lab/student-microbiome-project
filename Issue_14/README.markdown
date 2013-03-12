Beta diversity plots with a explicit time axis
==============================================

QIIME: `5275ab17da0f7858d66fc9c242f3335ff5b37881`
SMP: `369024113ea798b2ac4d9d9c23ed07fa680ff506`

Load the mapping file and unzip the distance matrices as we will need them.

```bash
load_remote_mapping_file.py -o StudentMicrobiomeProject-map.txt -k 0AvglGXLayhG7dDFUZ3JVVkFrTFFjMWJDWTZheVVROVE
BETADIV=../beta-diversity/
MAPF=StudentMicrobiomeProject-map.txt

# echo "gunzipping ..."
gunzip -c ${BETADIV}unweighted_unifrac_dm.txt.gz > ./unweighted_unifrac_dm.txt
gunzip -c ${BETADIV}weighted_unifrac_dm.txt.gz > ./weighted_unifrac_dm.txt
```

To extract the sample IDs that are viable for time-series analysis for each individual body site, use the following python code:

```python
from qiime.filter import sample_ids_from_metadata_description as s
mapping_file = open('StudentMicrobiomeProject-map.txt', 'U').readlines()
sids_tongue = s(mapping_file,'TongueTimeseries:Yes')
sids_forhead = s(mapping_file,'ForeheadTimeseries:Yes')
sids_palm = s(mapping_file,'PalmTimeseries:Yes')
sids_gut = s(mapping_file,'GutTimeseries:Yes')

fd = open('sample_ids_to_keep.txt', 'w')
# what we care about is the union of all these sets of sample identifiers
for sid in list(set(sids_tongue) | set(sids_forhead) | set(sids_palm) | set(sids_gut)):
    fd.write('%s\n' % sid)

fd.close()
```
Now that we have a list of sample identifiers we care about, let's filter the distance matrices and create the plots, note that it might take a while (10-15 mins) for the 3d plots to be created 8GB of memory should be enough.

```bash
echo "filtering the distance matrices ..."
filter_distance_matrix.py -i unweighted_unifrac_dm.txt -o unweighted_unifrac_dm.tsonly.txt --sample_id_fp=sample_ids_to_keep.txt -m ${MAPF}
filter_distance_matrix.py -i weighted_unifrac_dm.txt -o weighted_unifrac_dm.tsonly.txt --sample_id_fp=sample_ids_to_keep.txt -m ${MAPF}

echo "creating pc files ..."
principal_coordinates.py -i unweighted_unifrac_dm.tsonly.txt -o unweighted_unifrac_pc.txt
principal_coordinates.py -i weighted_unifrac_dm.tsonly.txt -o weighted_unifrac_pc.txt

echo "WeeksSinceStart vectors ..."
make_3d_plots.py -i unweighted_unifrac_pc.txt -m ${MAPF} -o all_sites_ts_only_unweighted_ordered_by_weekssincestart_vectors --custom_axes='WeeksSinceStart' --add_vectors='PersonalIDBodySite,WeeksSinceStart'
make_3d_plots.py -i weighted_unifrac_pc.txt -m ${MAPF} -o all_sites_ts_only_weighted_ordered_by_weekssincestart_vectors --custom_axes='WeeksSinceStart' --add_vectors='PersonalIDBodySite,WeeksSinceStart'
```




