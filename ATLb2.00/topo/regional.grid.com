#
set echo
#
# --- create a standard mercator regional.grid.[ab]
#
cd ~/hycom/ATLb2.00/topo
#
touch      fort.61
/bin/rm -f fort.61*
#
# --- temporary regional.grid.b
#
cat >! regional.grid.b <<'E-o-D'
  57	  'idm   ' = longitudinal array size
  52	  'jdm   ' = latitudinal  array size
'E-o-D'
setenv FOR061A fort.61A
../../ALL/topo/src/grid_mercator<<'E-o-D'
 57	'idm   ' = longitudinal array size
 52	'jdm   ' = latitudinal  array size
  0	'mapflg' = map flag (0=mercator,2=uniform,4=f-plane)
  1.0	'pntlon' = longitudinal reference grid point on pressure grid
263.0	'reflon' = longitude of reference grid point on pressure grid
  2.0	'grdlon' = longitudinal grid size (degrees)
 11.0	'pntlat' = latitudinal  reference grid point on pressure grid
  0.0	'reflat' = latitude of  reference grid point on pressure grid
  2.0	'grdlat' = latitudinal  grid size at the equator (degrees)
'E-o-D'
mv fort.61  regional.grid.b
mv fort.61A regional.grid.a
