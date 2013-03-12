#!/bin/bash
#Yoshiki Vazquez-Baeza
#QIIME 55c780f549fedacb8d4f0050b9ce99c73934cc49
#Mon Mar 11 18:24:03 MDT 2013

for PERSONAL_ID in `cat personal_ids.txt`
	do
	echo "The curent personal ID is ${PERSONAL_ID}"	
	echo "filtering ..." > job_${PERSONAL_ID}.txt
	echo "filter_distance_matrix.py -m ${PWD}/mapping_file.txt -s 'PersonalID:${PERSONAL_ID};AnyTimeseries:Yes' -o ${PWD}/per_subject_dm/${PERSONAL_ID}_unweighted_dm.txt -i ${PWD}/unweighted_unifrac_dm.txt" >> job_${PERSONAL_ID}.txt
	echo "principal coordinating ..." >> job_${PERSONAL_ID}.txt
	echo "principal_coordinates.py -i ${PWD}/per_subject_dm/${PERSONAL_ID}_unweighted_dm.txt -o ${PWD}/per_subject_pc/${PERSONAL_ID}_unweighted_pc.txt" >> job_${PERSONAL_ID}.txt
	echo "making 3d plots ..." >> job_${PERSONAL_ID}.txt
	echo "make_3d_plots.py -i ${PWD}/per_subject_pc/${PERSONAL_ID}_unweighted_pc.txt -m ${PWD}/mapping_file.txt -o ${PWD}/per_subject_3d_plots/${PERSONAL_ID}_unweighted --add_vectors='PersonalIDBodySite,WeeksSinceStart'" >> job_${PERSONAL_ID}.txt

	echo "filtering ..." >> job_${PERSONAL_ID}.txt
	echo "filter_distance_matrix.py -m ${PWD}/mapping_file.txt -s 'PersonalID:${PERSONAL_ID};AnyTimeseries:Yes' -o ${PWD}/per_subject_dm/${PERSONAL_ID}_weighted_dm.txt -i ${PWD}/weighted_unifrac_dm.txt " >> job_${PERSONAL_ID}.txt
	echo "principal coordinating ..." >> job_${PERSONAL_ID}.txt
	echo "principal_coordinates.py -i ${PWD}/per_subject_dm/${PERSONAL_ID}_weighted_dm.txt -o ${PWD}/per_subject_pc/${PERSONAL_ID}_weighted_pc.txt" >> job_${PERSONAL_ID}.txt
	echo "making 3d plots ..." >> job_${PERSONAL_ID}.txt
	echo "make_3d_plots.py -i ${PWD}/per_subject_pc/${PERSONAL_ID}_weighted_pc.txt -m ${PWD}/mapping_file.txt -o ${PWD}/per_subject_3d_plots/${PERSONAL_ID}_weighted --add_vectors='PersonalIDBodySite,WeeksSinceStart'" >> job_${PERSONAL_ID}.txt
	
	echo "sh job_${PERSONAL_ID}.txt" > boom.txt
	cluster_jobs_8.py -ms boom.txt ${PERSONAL_ID}

	done

	
