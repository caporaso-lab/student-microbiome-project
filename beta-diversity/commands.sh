#!/bin/bash
# QIIME 1.6.0
# SMP: c5161db845d9265d0d76d055daada75bc7ae4afe
# map.txt MD5: d8c3a8dc7246b0319c895967a068a485
# Author: Jai Ram Rideout

# Generate per-body-site PCoA plots (2D and 3D) using time-series-only data for
# weighted and unweighted UniFrac. Color by PersonalID.
MAP=map.txt
PID=PersonalID
UNIV=University
for ZIPFILE in *weighted_unifrac_dm.*_ts_only_pc.txt.gz
do
PCFILE=$(echo $ZIPFILE | sed 's/.gz//')
  OUTDIRBASE1=$(echo $PCFILE | sed 's/_pc.txt//')_personal_id_plots
  OUTDIRBASE2=$(echo $PCFILE | sed 's/_pc.txt//')_University_plots
  mkdir $OUTDIRBASE1
  mkdir $OUTDIRBASE2
  gunzip -c $ZIPFILE > $PCFILE
  make_2d_plots.py -i $PCFILE -m $MAP -b $UNIV -o $OUTDIRBASE1/2d_discrete_personal_id
  make_3d_plots.py -i $PCFILE -m $MAP -b $UNIV -o $OUTDIRBASE1/3d_discrete_personal_id -s scaled,unscaled
  make_2d_plots.py -i $PCFILE -m $MAP -b $UNIV -o $OUTDIRBASE2/2d_discrete_university
  make_3d_plots.py -i $PCFILE -m $MAP -b $UNIV -o $OUTDIRBASE2/3d_discrete_university -s scaled,unscaled
  tar czf $OUTDIRBASE1.tar.gz $OUTDIRBASE1
  tar czf $OUTDIRBASE2.tar.gz $OUTDIRBASE2
done