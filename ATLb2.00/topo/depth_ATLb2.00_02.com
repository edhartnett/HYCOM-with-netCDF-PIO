#
set echo
#
# --- create a once smoothed hycom topography.
#
setenv FOR051  fort.51
setenv FOR051A fort.51A
setenv FOR061  fort.61
setenv FOR061A fort.61A
#
/bin/rm -f $FOR051 $FOR051A $FOR061 $FOR061A
#
ln -s depth_ATLb2.00_01.b fort.51
ln -s depth_ATLb2.00_01.a fort.51A
#
../../ALL/topo/src/smooth<<'E-o-D'
1
'E-o-D'
mv fort.61  depth_ATLb2.00_02.b
mv fort.61A depth_ATLb2.00_02.a
#
/bin/rm -f $FOR051 $FOR051A
