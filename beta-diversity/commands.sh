#!/bin/bash
# QIIME:30d29ee4f6b5075e515373dbb08548a432315fb3
# SMP:104f2ab1515052775c47b200c839cf28f3d66573
# generate per body site principal coordinates files

rm input_dms.txt
echo "unweighted_unifrac_dm.forehead_ts_only.txt.gz" >> input_dms.txt
echo "unweighted_unifrac_dm.gut_ts_only.txt.gz" >> input_dms.txt
echo "unweighted_unifrac_dm.palm_ts_only.txt.gz" >> input_dms.txt
echo "unweighted_unifrac_dm.tongue_ts_only.txt.gz" >> input_dms.txt
echo "weighted_unifrac_dm.forehead_ts_only.txt.gz" >> input_dms.txt
echo "weighted_unifrac_dm.gut_ts_only.txt.gz" >> input_dms.txt
echo "weighted_unifrac_dm.palm_ts_only.txt.gz" >> input_dms.txt
echo "weighted_unifrac_dm.tongue_ts_only.txt.gz" >> input_dms.txt

for i in `cat input_dms.txt`
    do
        FILENAME=$(echo $i | sed 's/.gz//')
        PC_FILENAME=pcoa/$(echo "${FILENAME}" | sed 's/.txt/_pc.txt/')
        echo "processing ${FILENAME}"
        echo "result is ${PC_FILENAME}"
        #continue
        gunzip -c $i > ${FILENAME}
        principal_coordinates.py -i ${FILENAME} -o ${PC_FILENAME}
        gzip -c ${PC_FILENAME} > ${PC_FILENAME}.gz
    done