#
set echo
#
# --- interpolate 5-min bathymetry to hycom landsea mask.
#
setenv FOR051  fort.51
setenv FOR061A fort.61A
#
/bin/rm -f $FOR051 $FOR061A
#
/bin/ln -s ~/topo_ieee/ds759.2/tbase.a fort.51
#
../../ALL/topo/src/landsea_05min <<'E-o-D'
 &TOPOG
  INTERP = -5,
  MTYPE  =  0,
 /
'E-o-D'
#
/bin/mv fort.61A landsea_IASb0.50.a
#
/bin/rm -f fort.51
