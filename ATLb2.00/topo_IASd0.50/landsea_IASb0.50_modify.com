#
set echo
#
# --- modify a hycom land/sea mask.
#
cd ~/hycom/IASb0.50/topo
#
/bin/rm -f fort.[56]1*
#
setenv FOR051A fort.51A
setenv FOR061A fort.61A
#
/bin/mv landsea_IASb0.50.a landsea_IASb0.50.a_orig
/bin/ln -s landsea_IASb0.50.a_orig fort.51A
../../ALL/topo/src/mask_modify <<'E-o-D'
7
0  35  35  10  10
0  38  38  10  10
0  36  36  11  11
0  40  40  35  35
0  29  32  39  39
0  33  33  46  46
0  32  33  47  47
'E-o-D'
/bin/mv fort.61A landsea_IASb0.50.a
#
/bin/rm -f fort.[56]1*
