#
set echo
#
# --- landmask a hycom topography.
#
cd ~/hycom/IASb0.50/topo
#
/bin/rm -f fort.[56]1*
#
setenv FOR051  fort.51
setenv FOR051A fort.51A
setenv FOR061  fort.61
setenv FOR061A fort.61A
#
/bin/ln -s depth_IASb0.50_01.b fort.51
/bin/ln -s depth_IASb0.50_01.a fort.51A
../../ALL/topo/src/topo_landmask <<'E-o-D'
landmask the Pacific and Cuba
8
 1 40  1 10
 1 17 11 25
19 27 11 19
28 28 11 11
33 34 40 40
31 37 39 39
36 36 38 38
37 37 37 37
'E-o-D'
/bin/mv fort.61  depth_IASb0.50_02.b
/bin/mv fort.61A depth_IASb0.50_02.a
#
/bin/rm -f fort.[56]1*
