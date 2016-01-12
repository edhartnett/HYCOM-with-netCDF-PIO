#
set echo
#
# --- map a HYCOM bathymetry.
#
cd ~/hycom/IASb0.50/topo
#
setenv FOR051A fort.51A
/bin/rm -f fort.[5]1*
ln -s depth_IASb0.50_01.b  fort.51
ln -s depth_IASb0.50_01.a  fort.51A
if (-e landsea_IASb0.50.a) then
  /bin/rm -f regional.mask.a
  /bin/ln -s landsea_IASb0.50.a regional.mask.a
endif
../../ALL/topo/src/topo_map
