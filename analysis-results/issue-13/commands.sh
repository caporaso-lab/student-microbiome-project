#!/bin/bash
# QIIME 1.6.0
# SMP: c5161db845d9265d0d76d055daada75bc7ae4afe
# map.txt MD5: 38b71653eb3f88de56ba2785569f7ef0
# Author: Jai Ram Rideout

# Run several methods in compare_categories.py over per-body-site
# time-series-only distance matrices using the PersonalID category.
DMDIR=../../beta-diversity
MAP=map.txt
PID=PersonalID
NUMPERMS=999
for ZIPFILE in $DMDIR/*weighted_unifrac_dm.*_ts_only.txt.gz
do
  DMFILE=$(echo $ZIPFILE | sed 's/.gz//')
  OUTDIR=$(basename $DMFILE | sed 's/.txt//')_stats
  mkdir $OUTDIR
  gunzip -c $ZIPFILE > $DMFILE
  compare_categories.py -i $DMFILE -m $MAP -c $PID -o $OUTDIR -n $NUMPERMS --method adonis
  compare_categories.py -i $DMFILE -m $MAP -c $PID -o $OUTDIR -n $NUMPERMS --method anosim
  compare_categories.py -i $DMFILE -m $MAP -c $PID -o $OUTDIR -n $NUMPERMS --method permdisp
done
