#
set echo
#
# --- create the relaxation blkdat.input from the expt version.
#
source EXPT.src
#
echo "Levitus (NOAA World Ocean Atlas 1994) Climatology"        >! blkdat.input
echo "  00	  'month ' = month of climatology (01 to 12)"   >> blkdat.input
egrep "'iversn'|'iexpt '|'yrflag'" ../../expt_${X}/blkdat.input >> blkdat.input
egrep "'idm   '|'jdm   '"          ../../expt_${X}/blkdat.input >> blkdat.input
egrep "'itest '|'jtest '|'kdm   '" ../../expt_${X}/blkdat.input >> blkdat.input
egrep "'nhybrd'|'nsigma'"          ../../expt_${X}/blkdat.input >> blkdat.input
egrep "'dp00. '|'ds00. '"          ../../expt_${X}/blkdat.input >> blkdat.input
egrep "'thflag'|'thbase'|'sigma '" ../../expt_${X}/blkdat.input >> blkdat.input
egrep "'thkmin'"                   ../../expt_${X}/blkdat.input >> blkdat.input
