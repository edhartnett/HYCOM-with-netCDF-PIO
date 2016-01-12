#
set echo
#
# --- interpolate 5-min to hycom bathymetry.
#
setenv FOR051  fort.51
setenv FOR061  fort.61
setenv FOR061A fort.61A
#
/bin/rm -f $FOR051 $FOR061 $FOR061A
#
/bin/ln -s ~/topo_ieee/ds759.2/tbase.a fort.51
#
../../ALL/topo/src/bathy_05min <<'E-o-D'
 &TOPOG
  COAST  = 20.0,
  INTERP = -5,
  MTYPE  =  0,
 /
'E-o-D'
#
/bin/mv fort.61  depth_IASb0.50_01.b
/bin/mv fort.61A depth_IASb0.50_01.a
#
/bin/rm -f fort.51
