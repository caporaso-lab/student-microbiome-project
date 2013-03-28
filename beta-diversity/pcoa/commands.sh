#!/bin/bash
# QIIME 1.6.0
# SMP: c5161db845d9265d0d76d055daada75bc7ae4afe
# map.txt MD5: d8c3a8dc7246b0319c895967a068a485
# Author: Jai Ram Rideout

# Generate per-body-site PCoA plots (2D and 3D) using time-series-only data for
# weighted and unweighted UniFrac. Color by PersonalID.
MAP=map.txt
PID=PersonalID
for ZIPFILE in *weighted_unifrac_dm.*_ts_only_pc.txt.gz
do
  PCFILE=$(echo $ZIPFILE | sed 's/.gz//')
  OUTDIRBASE=$(echo $PCFILE | sed 's/_pc.txt//')_plots
  mkdir $OUTDIRBASE
  gunzip -c $ZIPFILE > $PCFILE
  make_2d_plots.py -i $PCFILE -m $MAP -b $PID -o $OUTDIRBASE/2d_discrete
  make_3d_plots.py -i $PCFILE -m $MAP -b $PID -o $OUTDIRBASE/3d_discrete -s scaled,unscaled
  tar czf $OUTDIRBASE.tar.gz $OUTDIRBASE
done
