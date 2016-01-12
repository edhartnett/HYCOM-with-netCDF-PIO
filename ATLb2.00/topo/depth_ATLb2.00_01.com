#
set echo
#
# --- convert micom to hycom topography.
# --- rotated coordinates.
#
setenv FOR051  fort.51
setenv FOR061  fort.61
setenv FOR061A fort.61A
#
/bin/rm -f $FOR051 $FOR061 $FOR061A
#
ln -s depth.51x56 fort.51
#
# --- two possible micom pakk formats (use topo_m2h or topo_malt2h).
#
../../ALL/topo/src/topo_m2h <<'E-o-D'
North Atlantic topography from ETOPO5 5-minute data set.
56 x 51 depth values (idm,jdm = 57,52).
2.0000 deg. mercator resolution, with equat=11.000.
lon: 97.000W to 13.000E.  lat: 19.606S to 62.195N.

'E-o-D'
mv fort.61  depth_ATLb2.00_01.b
mv fort.61A depth_ATLb2.00_01.a
#
/bin/rm -f $FOR051
