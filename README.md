student-microbiome-project
==========================

Central repository for data and analysis tools for the Student Microbiome Project (SMP). This repository will store data and code related to the Student Microbiome Project, and will be kept private until publication of the study. Related data include the [personal microbiome delivery system](https://github.com/qiime/personal-microbiome-delivery) being developed in Greg Caporaso's lab by John Chase, and the [Student Microbiome Project metadata mapping file](https://docs.google.com/spreadsheet/ccc?key=0AvglGXLayhG7dDFUZ3JVVkFrTFFjMWJDWTZheVVROVE). 

Description of data files
-------------------------

``otu_tables/`` : closed-reference OTU tables

``otu_table_stats/`` : output of running ``per_library_stats.py`` on all files in ``otu_tables/``

Analysis notes
--------------

Original mapping files had some overlapping sample IDs (where p, f, g, and t are palm, forehead, gut and tongue, respectively):

```
In [9]: set(p.keys()) & set(f.keys())
Out[9]: set(['ExtB3.1', 'ExtB5.1', 'ExtBS3.1', 'ExtBS3.2'])

In [10]: set(g.keys()) & set(t.keys())
Out[10]: set(['ExtB2.1', 'ExtBS2.2', 'ExtB4.1'])
```

These were renamed in the per-body-site mapping files to, for example, ``ExtB3.p.1`` and ``ExtB3.f.1`` to support a single mapping file to rule them all. 

After merging mapping files, the following steps were taken to generate information on the number of weeks of samples that were available for each subject.

```
from qiime.parse import parse_mapping_file_to_dict

s, _ = parse_mapping_file_to_dict(open('./smp_mapping_file_all_weeks.tsv','U'))

sid_to_ind = {}

ind_to_weeks = {}

for k,v in s.items():       ind = v['PersonalID']
    week = v['Week']       try:
        ind_to_weeks[ind].append(week)                                                                                                                                                  except KeyError:
        ind_to_weeks[ind] = [week]
    sid_to_ind[k] = ind

for k,v in sid_to_ind.items():
    if 'na' in ind_to_weeks[v]:
        print '\t'.join([k,'na','na','na'])
    else:
        print '\t'.join([k,str(len(set(ind_to_weeks[v]))),s[k]['Week'].replace('week.',''),str(min(map(float,[e.replace('week.','') for e in ind_to_weeks[v]]))), str(max(map(float,[e.replace('week.','') for e in ind_to_weeks[v]])))])
```

Demultiplexing
--------------
```
echo "split_libraries_fastq.py -i /home/shared/Illumina_hiseq_UCB/Illumina071712/Project_Fierer_71712/Sample_Fierer_16sV4_SMP1_tongue_71712/Fierer_16sV4_SMP1_tongue_71712_NoIndex_L001_R1_001.fastq.gz -b /home/shared/Illumina_hiseq_UCB/Illumina071712/Project_Fierer_71712/Sample_Fierer_16sV4_SMP1_tongue_71712/Fierer_16sV4_SMP1_tongue_71712_NoIndex_L001_R2_001.fastq.gz -m /home/caporaso/analysis/student_microbiome/tongue_all.tsv --rev_comp_mapping_barcodes -o /home/caporaso/analysis/student_microbiome/slout_tongue/" | qsub -keo -N stu_tongue

echo "split_libraries_fastq.py -i /home/shared/Illumina_hiseq_UCB/Illumina071712/Project_Fierer_71712/Sample_Fierer_16sV4_SMP2_palm_71712/Fierer_16sV4_SMP2_palm_71712_NoIndex_L002_R1_001.fastq.gz -b /home/shared/Illumina_hiseq_UCB/Illumina071712/Project_Fierer_71712/Sample_Fierer_16sV4_SMP2_palm_71712/Fierer_16sV4_SMP2_palm_71712_NoIndex_L002_R2_001.fastq.gz -m /home/caporaso/analysis/student_microbiome/palm_all.tsv --rev_comp_mapping_barcodes -o /home/caporaso/analysis/student_microbiome/slout_palm/" | qsub -keo -N stu_palm

echo "split_libraries_fastq.py -i /home/shared/Illumina_hiseq_UCB/Illumina071712/Project_Fierer_71712/Sample_Fierer_16sV4_SMP3_forehead_71712/Fierer_16sV4_SMP3_forehead_71712_NoIndex_L003_R1_001.fastq.gz -b /home/shared/Illumina_hiseq_UCB/Illumina071712/Project_Fierer_71712/Sample_Fierer_16sV4_SMP3_forehead_71712/Fierer_16sV4_SMP3_forehead_71712_NoIndex_L003_R2_001.fastq.gz -m /home/caporaso/analysis/student_microbiome/forehead_all.tsv --rev_comp_mapping_barcodes -o /home/caporaso/analysis/student_microbiome/slout_forehead/" | qsub -keo -N stu_forehead

echo "split_libraries_fastq.py -i /home/shared/Illumina_hiseq_UCB/Illumina071712/Project_Fierer_71712/Sample_Fierer_16sV4_SMP4_fecal_71712/Fierer_16sV4_SMP4_fecal_71712_NoIndex_L004_R1_001.fastq.gz -b /home/shared/Illumina_hiseq_UCB/Illumina071712/Project_Fierer_71712/Sample_Fierer_16sV4_SMP4_fecal_71712/Fierer_16sV4_SMP4_fecal_71712_NoIndex_L004_R2_001.fastq.gz -m /home/caporaso/analysis/student_microbiome/gut_all.tsv --rev_comp_mapping_barcodes -o /home/caporaso/analysis/student_microbiome/slout_fecal/" | qsub -keo -N stu_fecal

# week one only data
echo "split_libraries_fastq.py -i /home/shared/Illumina_hiseq_UCB/Illumina03022012/Project_Fierer_030212/Sample_NF_SM_H/NF_SM_H_NoIndex_L001_R1_001.fastq -b /home/shared/Illumina_hiseq_UCB/Illumina03022012/Project_Fierer_030212/Sample_NF_SM_H/NF_SM_H_NoIndex_L001_R2_001.fastq -m /scratch/caporaso/student_microbiome/wk1only_map.txt -o /scratch/caporaso/student_microbiome/slout_wk1_only/ --rev_comp_mapping_barcodes" | qsub -keo -N student_sl
```

OTU picking
-----------
```
echo "pick_reference_otus_through_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/slout_fecal/seqs.fna -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_fecal/ -r /Users/caporaso/data/gg_12_10_otus/rep_set/97_otus.fasta -t /Users/caporaso/data/gg_12_10_otus/taxonomy/97_otu_taxonomy.txt -aO 50 -p /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_params.txt" | qsub -keo -N smp-fecal -l pvmem=16gb -q memroute

echo "pick_reference_otus_through_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/slout_forehead/seqs.fna -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_forehead/ -r /Users/caporaso/data/gg_12_10_otus/rep_set/97_otus.fasta -t /Users/caporaso/data/gg_12_10_otus/taxonomy/97_otu_taxonomy.txt -aO 50 -p /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_params.txt" | qsub -keo -N smp-forehead -l pvmem=16gb -q memroute

echo "pick_reference_otus_through_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/slout_tongue/seqs.fna -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_tongue/ -r /Users/caporaso/data/gg_12_10_otus/rep_set/97_otus.fasta -t /Users/caporaso/data/gg_12_10_otus/taxonomy/97_otu_taxonomy.txt -aO 50 -p /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_params.txt" | qsub -keo -N smp-tongue -l pvmem=16gb -q memroute

echo "pick_reference_otus_through_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/slout_palm/seqs.fna -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_palm/ -r /Users/caporaso/data/gg_12_10_otus/rep_set/97_otus.fasta -t /Users/caporaso/data/gg_12_10_otus/taxonomy/97_otu_taxonomy.txt -aO 50 -p /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_params.txt" | qsub -keo -N smp-palm -l pvmem=16gb -q memroute

# week one only data
echo "pick_reference_otus_through_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/slout_wk1_only/seqs.fna -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_wk1/ -r /Users/caporaso/data/gg_12_10_otus/rep_set/97_otus.fasta -t /Users/caporaso/data/gg_12_10_otus/taxonomy/97_otu_taxonomy.txt -aO 50 -p /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_params.txt" | qsub -keo -N smp-wk1 -l pvmem=16gb -q memroute

```


OTU table filtering
-------------------
First we remove samples that are not in the mapping file. Some of these had issues during prep, etc.
```
filter_samples_from_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_fecal/uclust_ref_picked_otus/otu_table.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_fecal/uclust_ref_picked_otus/otu_table_sfiltered.biom --sample_id_fp /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv 

filter_samples_from_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_forehead/uclust_ref_picked_otus/otu_table.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_forehead/uclust_ref_picked_otus/otu_table_sfiltered.biom --sample_id_fp /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv 

filter_samples_from_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_tongue/uclust_ref_picked_otus/otu_table.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_tongue/uclust_ref_picked_otus/otu_table_sfiltered.biom --sample_id_fp /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv 

filter_samples_from_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_palm/uclust_ref_picked_otus/otu_table.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_palm/uclust_ref_picked_otus/otu_table_sfiltered.biom --sample_id_fp /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv
```

Next we create a master OTU table by merging the four body sites
```
echo "merge_otu_tables.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_fecal/uclust_ref_picked_otus/otu_table_sfiltered.biom,/scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_forehead/uclust_ref_picked_otus/otu_table_sfiltered.biom,/scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_tongue/uclust_ref_picked_otus/otu_table_sfiltered.biom,/scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC_palm/uclust_ref_picked_otus/otu_table_sfiltered.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/all_body_sites_sfiltered.biom" | qsub -keo -N merge-smp -l pvmem=32gb -q memroute
```

Another round of filtering. We're requring a minimum of 10k sequences per sample.
```
echo "filter_samples_from_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/all_body_sites_sfiltered.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/all_body_sites_sfiltered_mc10k.biom -n 10000" | qsub -keo -N smp-filt -l pvmem=32gb -q memroute
```

Next we apply a filtering step to remove OTUs that show up in high abundance in our negative control samples (i.e., extraction and swab blanks). The specific revision of the script used in this step can be found [here](https://gist.github.com/4197809/bf1df4cae8f742777736ee08d3c1c4471753996a). Then filter samples which have <10k sequences after this OTU filtering step, and finally remove the negative control samples from the OTU table.
```
echo 'cd /Users/caporaso/code/filter_control_otus ; python filter_control_otus.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/all_body_sites_sfiltered_mc10k.biom -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -a 0.005 -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/control_filtering/ -s "PersonalID:SwabBlank,NTC,ntc,ExtBlank"' | qsub -keo -N smp-filt -l pvmem=64gb -q memroute

echo "filter_samples_from_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/control_filtering/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/control_filtering/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k.biom -n 10000" | qsub -keo -N smp-filt -l pvmem=32gb -q memroute

echo 'filter_samples_from_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/control_filtering/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/control_filtering/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered.biom -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -s "PersonalID:*,!SwabBlank,!NTC,!ntc,!ExtBlank"' | qsub -keo -N smp-filt -l pvmem=32gb -q memroute
```

Now identify samples that are likely to be mislabeled and remove them.
```
# filter OTUs in fewer than 10 samples
echo "filter_otus_from_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/control_filtering/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10.biom -s 10" | qsub -keo -N filtotus -l pvmem=32gb -q memroute

# rarify to 1000 sequences/sample
echo "single_rarefaction.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10_even1000.biom -d 1000" | qsub -keo -N filtotus -l pvmem=32gb -q memroute

# identify mislabeled samples
echo "supervised_learning.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10_even1000.biom -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -c "BodySite" -e cv5 -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/s10_1000/" | qsub -keo -N smp_rf  -l pvmem=16gb -q memroute

# remove the mislabeled samples
echo "filter_samples_from_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/control_filtering/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_no_mislabeled_010.biom -m /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/s10_1000/mislabeling.txt -s 'mislabeled_probability_above_0.90:FALSE'" | qsub -keo -N filtmisl -l pvmem=32gb -q memroute

# the output of the previous step becomes our master OTU table for diversity analyses
cp /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_no_mislabeled_010.biom /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table.biom

# SMP mapping fill pulled from GDocs at ~8:55pm MT 12/20/12

per_library_stats.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table.biom -o /scratch/caporaso/student_microbiome/SMP-map.out.tsv -m /scratch/caporaso/student_microbiome/SMP-map.tsv

```

Diversity analyses
------------------

Compute UniFrac distance matrices and generate PCoA plots.
```
echo "single_rarefaction.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table_even1000.biom -d 10000" | qsub -keo -N smp-srare

# Note: the output file was incorrect in the above command (should have made it otu_table_even10000.biom, so 10k instead of 1k). The 
# error follows through in the next couple of commands.
parallel_beta_diversity.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table_even1000.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/ -O 100 -U /Users/caporaso/bin/cluster_jobs_8.py -t /Users/caporaso/data/gg_12_10_otus/trees/97_otus.tree

cd /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/
mv unweighted_unifrac_otu_table_even1000.txt unweighted_unifrac_dm.txt
mv weighted_unifrac_otu_table_even1000.txt weighted_unifrac_dm.txt

parallel_multiple_rarefactions.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table.biom -m 10 -x 10000 -s 999 -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/arare_10000//rarefaction/  --jobs_to_start 50

parallel_alpha_diversity.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/arare_10000//rarefaction/ -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/arare_10000//alpha_div/ --metrics PD_whole_tree,chao1,observed_species,shannon -t /Users/caporaso/data/gg_12_10_otus/trees/97_otus.tree --jobs_to_start 50

collate_alpha.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/arare_10000/alpha_div -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/arare_10000/alpha_div_collated/

```

Generate timeseries only full and per-BodySite OTU tables
---------------------------------------------------------
```

filter_samples_from_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table.ts_only.biom -s "AnyTimeseries:Yes"  -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv

split_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table.ts_only.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table_ts_by_body_site/ -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv -f BodySite

filter_samples_from_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project/otu-tables/otu_table_even10000.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project/otu-tables/otu_table_even10000.ts_only.biom -s "AnyTimeseries:Yes"  -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv

split_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project/otu-tables/otu_table_even10000.ts_only.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table_even_10000_ts_by_body_site/ -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv -f BodySite
```

Filter distance matrices to timeseries only
-------------------------------------------
```
echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/unweighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/unweighted_unifrac_dm.ts_only.txt -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv -s "AnyTimeseries:Yes"' | qsub -keo -N dmts -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/unweighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/unweighted_unifrac_dm.gut_ts_only.txt -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv -s "GutTimeseries:Yes"' | qsub -keo -N dmgutts -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/unweighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/unweighted_unifrac_dm.palm_ts_only.txt -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv -s "PalmTimeseries:Yes"' | qsub -keo -N dmpalmts -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/unweighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/unweighted_unifrac_dm.forehead_ts_only.txt -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv -s "ForeheadTimeseries:Yes"' | qsub -keo -N dmforets -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/unweighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/unweighted_unifrac_dm.tongue_ts_only.txt -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv -s "TongueTimeseries:Yes"' | qsub -keo -N dmtonguets -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/weighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/weighted_unifrac_dm.ts_only.txt -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv -s "AnyTimeseries:Yes"' | qsub -keo -N dmts -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/weighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/weighted_unifrac_dm.gut_ts_only.txt -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv -s "GutTimeseries:Yes"' | qsub -keo -N dmgutts -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/weighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/weighted_unifrac_dm.palm_ts_only.txt -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv -s "PalmTimeseries:Yes"' | qsub -keo -N dmpalmts -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/weighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/weighted_unifrac_dm.forehead_ts_only.txt -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv -s "ForeheadTimeseries:Yes"' | qsub -keo -N dmforets -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/weighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/weighted_unifrac_dm.tongue_ts_only.txt -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv -s "TongueTimeseries:Yes"' | qsub -keo -N dmtonguets -l pvmem=8gb -q memroute
```

Run PCoA on select distance matrices
------------------------------------

```
echo 'principal_coordinates.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/unweighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/unweighted_unifrac_pc.txt' | qsub -keo -N pc -l pvmem=8gb -q memroute

echo 'principal_coordinates.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/unweighted_unifrac_dm.ts_only.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/unweighted_unifrac_pc.ts_only.txt' | qsub -keo -N pcts -l pvmem=8gb -q memroute

echo 'principal_coordinates.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/weighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/weighted_unifrac_pc.txt' | qsub -keo -N pc -l pvmem=8gb -q memroute

echo 'principal_coordinates.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/weighted_unifrac_dm.ts_only.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000/weighted_unifrac_pc.ts_only.txt' | qsub -keo -N pcts -l pvmem=8gb -q memroute


```

Generate alpha rarefaction data and plots
-----------------------------------------
```
echo "alpha_rarefaction.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/arare_10000/ -e 10000 -m /scratch/caporaso/student_microbiome/SMP-map_w_ts.tsv -p /scratch/caporaso/student_microbiome/student-microbiome-project/parameters/arare_params.txt -t /Users/caporaso/data/gg_12_10_otus/trees/97_otus.tree -aO 50" | qsub -k oe -N smp-arare -l pvmem=8gb -q memroute

```

