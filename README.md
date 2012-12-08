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

echo "single_rarefaction.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10_even1000.biom -d 1000" | qsub -keo -N filtotus -l pvmem=32gb -q memroute

echo "supervised_learning.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10_even1000.biom -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -c "BodySite" -e cv5 -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/s10_1000/" | qsub -keo -N smp_rf  -l pvmem=16gb -q memroute

echo "filter_samples_from_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10_no_mislabeled_010.biom -m /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/s10_1000/mislabeling.txt -s 'mislabeled_at_0.10:FALSE'" | qsub -keo -N filtmisl -l pvmem=32gb -q memroute
```

Diversity analyses
------------------

Compute UniFrac distance matrices and generate PCoA plots.
```
echo "beta_diversity_through_plots.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10_no_mislabeled_010.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/bdiv_even10000 -e 10000 -aO 100 -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -p /scratch/caporaso/student_microbiome/student-microbiome-project/parameters/bdiv_params.txt -t /Users/caporaso/data/gg_12_10_otus/trees/97_otus.tree" | qsub -keo -N smp-bdiv -l pvmem=32gb -q memroute
```

Create per-body-site distance matrices
```
echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/unweighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/unweighted_unifrac_forehead_dm.txt -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -s "BodySite:forehead"' | qsub -keo -N dmforehead -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/unweighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/unweighted_unifrac_gut_dm.txt -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -s "BodySite:gut"' | qsub -keo -N dmgut -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/unweighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/unweighted_unifrac_palm_dm.txt -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -s "BodySite:palm"' | qsub -keo -N dmpalm -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/unweighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/unweighted_unifrac_tongue_dm.txt -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -s "BodySite:tongue"' | qsub -keo -N dmtongue -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/weighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/weighted_unifrac_forehead_dm.txt -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -s "BodySite:forehead"' | qsub -keo -N dmforehead -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/weighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/weighted_unifrac_gut_dm.txt -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -s "BodySite:gut"' | qsub -keo -N dmgut -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/weighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/weighted_unifrac_palm_dm.txt -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -s "BodySite:palm"' | qsub -keo -N dmpalm -l pvmem=8gb -q memroute

echo 'filter_distance_matrix.py -i /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/weighted_unifrac_dm.txt -o /scratch/caporaso/student_microbiome/student-microbiome-project/beta-diversity/weighted_unifrac_tongue_dm.txt -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -s "BodySite:tongue"' | qsub -keo -N dmtongue -l pvmem=8gb -q memroute
```


Compute alpha diversity and rarefaction plots.
```
echo "alpha_rarefaction.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10_no_mislabeled_010.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/arare_even10000 -e 10000 -aO 100 -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -p /scratch/caporaso/student_microbiome/student-microbiome-project/parameters/bdiv_params.txt -t /Users/caporaso/data/gg_12_10_otus/trees/97_otus.tree" | qsub -keo -N smp-arare -l pvmem=32gb -q memroute

# Oops, forgot to include shannon in the set of metrics to apply so computing after.

parallel_alpha_diversity.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/arare_even10000/rarefaction/ -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/arare_even10000/alpha_div_s/ --jobs_to_start 100 -m shannon
collate_alpha.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/arare_even10000/alpha_div_s/ -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/arare_even10000/alpha_div_collated/


```

Split the full OTU table into per-body site OTU tables
```
echo "split_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10_no_mislabeled_010.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/ -m /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv -f BodySite" | qsub -keo -N split-smp -l pvmem=64gb -q memroute
```

Sort the OTU tables and build taxa summary plots. The sorted OTU tables are the final working per-body-site OTU tables.
```
echo "sort_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10_no_mislabeled_010_gut.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/otu_table_gut.biom -l /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv ; summarize_taxa_through_plots.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/otu_table_gut.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/taxa_plots_gut/" | qsub -keo -N gut-tax -l pvmem=32gb -q memroute

echo "sort_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10_no_mislabeled_010_tongue.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/otu_table_tongue.biom -l /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv ; summarize_taxa_through_plots.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/otu_table_tongue.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/taxa_plots_tongue/" | qsub -keo -N tongue-tax -l pvmem=32gb -q memroute

echo "sort_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10_no_mislabeled_010_forehead.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/otu_table_forehead.biom -l /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv ; summarize_taxa_through_plots.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/otu_table_forehead.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/taxa_plots_forehead/" | qsub -keo -N forehead-tax -l pvmem=32gb -q memroute

echo "sort_otu_table.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/all_body_sites_sfiltered_mc10k.control_filtered_median_0.005_mc10k_sfiltered_s10_no_mislabeled_010_palm.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/otu_table_palm.biom -l /scratch/caporaso/student_microbiome/StudentMicrobiomeProject-map.tsv ; summarize_taxa_through_plots.py -i /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/otu_table_palm.biom -o /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/per_BodySite_otu_tables/taxa_plots_palm/"  | qsub -keo -N palm-tax -l pvmem=32gb -q memroute

```
