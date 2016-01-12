#
set echo
#
# --- form subregion bathymetry file, ATLb2.00 to IASb0.50.
#
# --- ALL/bin/hycom_ij2lonlat can be used to find co-located points
# ---  on the two grids.  Since subregion p(1,1) must be on the
# ---  original grid, this is usually the point to reference.
# --- In this case:
# ---   hycom_ij2lonlat 1 1 ~/hycom/IASb0.50/topo/regional.grid.a
# ---    97.000W   3.997N
# ---   hycom_ij2lonlat 1 13 ~/hycom/ATLb2.00/topo/regional.grid.a
# ---    97.000W   3.997N
#
setenv D ../topo
setenv R ../../IASb0.50/topo
#
touch regional.grid.a
if (-z regional.grid.a) then
  ln -s $D/regional.grid.a .
  ln -s $D/regional.grid.b .
endif
#
touch   ${R}/depth_IASb0.50_99.[ab]
/bin/rm ${R}/depth_IASb0.50_99.[ab]
../../ALL/subregion/src/isub_topog <<E-o-D
${D}/depth_ATLb2.00_01.b
${R}/depth_IASb0.50_99.b
depth_ATLb2.00_01 subregioned to IASb0.50 via isub_topog
  93	  'idm   ' = longitudinal array size
  65	  'jdm   ' = latitudinal  array size
   1	  'irefi ' = longitudinal input  reference location
  13	  'jrefi ' = latitudinal  input  reference location
   1	  'irefo ' = longitudinal output reference location
   1	  'jrefo ' = latitudinal  output reference location
   4	  'ijgrd ' = integer scale factor between input and output grids
E-o-D
