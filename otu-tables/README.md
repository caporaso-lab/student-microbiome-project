OTU table
6414a681f5b128430deb03e2a964aa4a  otu_table.biom.gz

OTU table rarefied to 10000 seqs/sample
854900c42a4e9c1f3b3dcb681602c8a6  otu_table_even10000.biom.gz

OTU table only including samples that are included in a timeseries
8e52fb54ff95615da5862f7eb26f7c42  otu_table.ts_only.biom.gz

OTU table only including samples that are included in forehead timeseries
20a34382c081683362c1e39632b38c41  otu_table_forehead.ts_only.biom.gz

OTU table only including samples that are included in gut timeseries
e1443d9b31b3bc20b372c5addd47e9ff  otu_table_gut.ts_only.biom.gz

OTU table only including samples that are included in palm timeseries
32e807c57c0c3b2e49aa0c217ee5b539  otu_table_palm.ts_only.biom.gz

OTU table only including samples that are included in tongue timeseries
2f05ccd09902b7ad306bfa5959fd6b9e  otu_table_tongue.ts_only.biom.gz

Timeseries files
1d5c425f5aafa98a4d51f2f6494c653c  otu_table_even10000_forehead.ts_only.biom.gz
7425b9aae0a3e05f299c2635c0b697a1  otu_table_even10000_gut.ts_only.biom.gz
015724757aa7b262897fa007b28ab473  otu_table_even10000_palm.ts_only.biom.gz
9e05418e5848be4ae4a96484cf7e31db  otu_table_even10000_tongue.ts_only.biom.gz
436393f1a872c2afcfcad27d29476497  otu_table_even10000.ts_only.biom.gz


File tracking
-------------

These commands were all run from the directory containing this README.md.
```
gzip -c /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table.biom > otu_table.biom.gz
gzip -c /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table.ts_only.biom > otu_table.ts_only.biom.gz
gzip -c /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table_ts_by_body_site/otu_table.ts_only_palm.biom > otu_table_palm.ts_only.biom.gz
gzip -c /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table_ts_by_body_site/otu_table.ts_only_tongue.biom > otu_table_tongue.ts_only.biom.gz
gzip -c /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table_ts_by_body_site/otu_table.ts_only_gut.biom > otu_table_gut.ts_only.biom.gz
gzip -c /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table_ts_by_body_site/otu_table.ts_only_forehead.biom > otu_table_forehead.ts_only.biom.gz
# note that the initial file name was bad on this one:
gzip -c /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table_even1000.biom > otu_table_even10000.biom.gz
gzip otu_table_even10000.ts_only.biom
gzip -c /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table_even_10000_ts_by_body_site/otu_table_even10000.ts_only_gut.biom > otu_table_even10000_gut.ts_only.biom.gz
gzip -c /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table_even_10000_ts_by_body_site/otu_table_even10000.ts_only_palm.biom > otu_table_even10000_palm.ts_only.biom.gz
gzip -c /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table_even_10000_ts_by_body_site/otu_table_even10000.ts_only_tongue.biom > otu_table_even10000_tongue.ts_only.biom.gz
gzip -c /scratch/caporaso/student_microbiome/student-microbiome-project-raw-data/ucrC/mislabeling/otu_table_even_10000_ts_by_body_site/otu_table_even10000.ts_only_forehead.biom > otu_table_even10000_forehead.ts_only.biom.gz
```
