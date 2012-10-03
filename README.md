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

Beta diversity
--------------

```
echo "beta_diversity_through_plots.py -i /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/palm_closed_ref_filtered_otu_table.biom -o /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/palm_bdiv_even10000/ -t /scratch/caporaso/gg_otus_4feb2011/trees/gg_97_otus_4feb2011.tre -m /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/StudentMicrobiomeProject-map.tsv -e 10000 -p /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/tongue_bdiv_even10000/bdiv_params.txt" | qsub -keo -N smppalm10000 -l pvmem=8gb -q memroute

echo "beta_diversity_through_plots.py -i /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/tongue_closed_ref_filtered_otu_table.biom -o /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/tongue_bdiv_even10000/ -t /scratch/caporaso/gg_otus_4feb2011/trees/gg_97_otus_4feb2011.tre -m /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/StudentMicrobiomeProject-map.tsv -e 10000 -p /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/tongue_bdiv_even10000/bdiv_params.txt" | qsub -keo -N smptongue10000 -l pvmem=8gb -q memroute

echo "beta_diversity_through_plots.py -i /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/gut_closed_ref_filtered_otu_table.biom -o /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/gut_bdiv_even10000/ -t /scratch/caporaso/gg_otus_4feb2011/trees/gg_97_otus_4feb2011.tre -m /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/StudentMicrobiomeProject-map.tsv -e 10000 -p /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/tongue_bdiv_even10000/bdiv_params.txt" | qsub -keo -N smpgut10000 -l pvmem=8gb -q memroute

echo "beta_diversity_through_plots.py -i /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/forehead_closed_ref_filtered_otu_table.biom -o /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/forehead_bdiv_even10000/ -t /scratch/caporaso/gg_otus_4feb2011/trees/gg_97_otus_4feb2011.tre -m /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/StudentMicrobiomeProject-map.tsv -e 10000 -p /Users/caporaso/analysis/student_microbiome/24sept2012/student-microbiome-project/otu_tables/tongue_bdiv_even10000/bdiv_params.txt" | qsub -keo -N smpforehead10000 -l pvmem=8gb -q memroute
```
