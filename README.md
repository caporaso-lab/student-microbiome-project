student-microbiome-project
==========================

Central repository for data and analysis tools for the Student Microbiome Project (SMP). This repository will store data and code related to the Student Microbiome Project, and will be kept private until publication of the study. Related data include the [personal microbiome delivery system](https://github.com/qiime/personal-microbiome-delivery) being developed in Greg Caporaso's lab by John Chase, and the [Student Microbiome Project metadata mapping file](https://docs.google.com/spreadsheet/ccc?key=0AvglGXLayhG7dDFUZ3JVVkFrTFFjMWJDWTZheVVROVE). 

Description of data files
-------------------------

``otu_tables/`` : closed-reference OTU tables - filtered tables have been manually processed by Gilberto Flores to remove samples with fewer than 10,000 sequences and OTUs that were present in negative controls. 

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
pick_subsampled_reference_otus_through_otu_table.py -i /scratch/caporaso/student_microbiome/slout_wk1_only/seqs.fna,/scratch/caporaso/student_microbiome/slout_fecal/seqs.fna,/scratch/caporaso/student_microbiome/slout_tongue/seqs.fna,/scratch/caporaso/student_microbiome/slout_forehead/seqs.fna,/scratch/caporaso/student_microbiome/slout_palm/seqs.fna -r /scratch/caporaso/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta -o /scratch/caporaso/student_microbiome/ucrss_fast/ -p /scratch/caporaso/student_microbiome/ucrss_params.txt -n student.microbiome -aO 50
```

OTU table filtering
-------------------

NEED GILBERT'S NOTES

Alpha rarefaction
-----------------

```
echo "alpha_rarefaction.py -i /scratch/gifl2111/SMP/Working_OTU_Tables/Forehead/forehead_closed_ref_working_otu_table.biom -o /home/caporaso/analysis/student_microbiome/24sept2012/forehead/arare_max10000/ -t /scratch/caporaso/gg_otus_4feb2011/trees/gg_97_otus_4feb2011.tre -m /scratch/gifl2111/SMP/Working_OTU_Tables/Forehead/forehead_working_mf.txt -aO 25 -e 10000" | qsub -keo -N smpmax10000 -l pvmem=8gb -q memroute

echo "alpha_rarefaction.py -i /scratch/gifl2111/SMP/Working_OTU_Tables/Gut/gut_closed_ref_working_otu_table.biom -o /home/caporaso/analysis/student_microbiome/24sept2012/gut/arare_max10000/ -t /scratch/caporaso/gg_otus_4feb2011/trees/gg_97_otus_4feb2011.tre -m /scratch/gifl2111/SMP/Working_OTU_Tables/Gut/gut_working_mf.txt -aO 25 -e 10000" | qsub -keo -N smpmax10000 -l pvmem=8gb -q memroute

echo "alpha_rarefaction.py -i /scratch/gifl2111/SMP/Working_OTU_Tables/Tongue/tongue_closed_ref_working_otu_table.biom -o /home/caporaso/analysis/student_microbiome/24sept2012/tongue/arare_max10000/ -t /scratch/caporaso/gg_otus_4feb2011/trees/gg_97_otus_4feb2011.tre -m /scratch/gifl2111/SMP/Working_OTU_Tables/Tongue/tongue_working_mf.txt -aO 25 -e 10000" | qsub -keo -N smpmax10000 -l pvmem=8gb -q memroute

echo "alpha_rarefaction.py -i /scratch/gifl2111/SMP/Working_OTU_Tables/Palm/palm_closed_ref_working_otu_table.biom -o /home/caporaso/analysis/student_microbiome/24sept2012/palm/arare_max10000/ -t /scratch/caporaso/gg_otus_4feb2011/trees/gg_97_otus_4feb2011.tre -m /scratch/gifl2111/SMP/Working_OTU_Tables/Palm/palm_working_mf.txt -aO 25 -e 10000" | qsub -keo -N smpmax10000 -l pvmem=8gb -q memroute
```

Beta diversity
--------------

```
echo "beta_diversity_through_plots.py -i /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/palm_closed_ref_filtered_otu_table.biom -o /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/palm_bdiv_even10000/ -t /scratch/caporaso/gg_otus_4feb2011/trees/gg_97_otus_4feb2011.tre -m /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/StudentMicrobiomeProject-map.tsv -e 10000 -p /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/bdiv_params.txt" | qsub -keo -N smppalm10000 -l pvmem=8gb -q memroute

echo "beta_diversity_through_plots.py -i /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/tongue_closed_ref_filtered_otu_table.biom -o /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/tongue_bdiv_even10000/ -t /scratch/caporaso/gg_otus_4feb2011/trees/gg_97_otus_4feb2011.tre -m /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/StudentMicrobiomeProject-map.tsv -e 10000 -p /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/bdiv_params.txt" | qsub -keo -N smptongue10000 -l pvmem=8gb -q memroute

echo "beta_diversity_through_plots.py -i /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/gut_closed_ref_filtered_otu_table.biom -o /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/gut_bdiv_even10000/ -t /scratch/caporaso/gg_otus_4feb2011/trees/gg_97_otus_4feb2011.tre -m /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/StudentMicrobiomeProject-map.tsv -e 10000 -p /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/bdiv_params.txt" | qsub -keo -N smpgut10000 -l pvmem=8gb -q memroute

echo "beta_diversity_through_plots.py -i /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/forehead_closed_ref_filtered_otu_table.biom -o /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/forehead_bdiv_even10000/ -t /scratch/caporaso/gg_otus_4feb2011/trees/gg_97_otus_4feb2011.tre -m /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/StudentMicrobiomeProject-map.tsv -e 10000 -p /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/bdiv_params.txt" | qsub -keo -N smpforehead10000 -l pvmem=8gb -q memroute
```

Merge OTU tables for combined analysis
--------------------------------------

```
echo "merge_otu_tables.py -i /Users/caporaso/analysis/student-microbiome-project/otu_tables/forehead_closed_ref_otu_table.biom,/Users/caporaso/analysis/student-microbiome-project/otu_tables/gut_closed_ref_otu_table.biom,/Users/caporaso/analysis/student-microbiome-project/otu_tables/palm_closed_ref_otu_table.biom,/Users/caporaso/analysis/student-microbiome-project/otu_tables/tongue_closed_ref_otu_table.biom -o /Users/caporaso/analysis/student-microbiome-project/otu_tables/closed_ref_otu_table.biom" | qsub -keo -N mergesmp -l pvmem=64gb -q memroute

echo "per_library_stats.py -i /Users/caporaso/analysis/student-microbiome-project/otu_tables/closed_ref_otu_table.biom > /Users/caporaso/analysis/student-microbiome-project/otu_table_stats/closed_ref_otu_table_per_lib_stats.txt" | qsub -keo -N smp_pls -l pvmem=32gb -q memroute

echo "beta_diversity_through_plots.py -i /Users/caporaso/analysis/student-microbiome-project/otu_tables/closed_ref_otu_table.biom -o /Users/caporaso/analysis/student-microbiome-project/bdiv_even10000/ -p /Users/caporaso/analysis/student-microbiome-project/parameters/bdiv_params.txt -t /Users/caporaso/data/gg_otus_4feb2011/trees/gg_97_otus_4feb2011.tre -m /Users/caporaso/analysis/student-microbiome-project/StudentMicrobiomeProject-map.tsv -e 10000 -aO 50" | qsub -keo -N smp_b10000 -l pvmem=8gb -q memroute
```
