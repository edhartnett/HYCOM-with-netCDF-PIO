#
set echo
#
# --- create the zonal clim blkdat.input from any expt version.
# --- unchanging part of input is from blkdat.template.
#
echo "Levitus (NOAA World Ocean Atlas 1994) Climatology"        >! blkdat.input
echo "  00	  'month ' = month of climatology (01 to 12)"   >> blkdat.input
egrep "'iversn'|'iexpt '|'yrflag'" ../../expt_02.1/blkdat.input >> blkdat.input
egrep "'mapflg'|'idm   '|'jdm   '" ../../expt_02.1/blkdat.input >> blkdat.input
egrep "'pntlat'|'reflat'|'grdlat'" ../../expt_02.1/blkdat.input >> blkdat.input
                                            cat blkdat.template >> blkdat.input
