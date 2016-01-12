#
set echo
#
# --- generate 2-d equal ocean processor partitions
#
cd ~/hycom/ATLb2.00/topo/partit
#
ln -s ../regional.grid.a .
ln -s ../regional.grid.b .
#
setenv FOR051  ../depth_ATLb2.00_01.b
setenv FOR051A ../depth_ATLb2.00_01.a
#
# --- 002 (1d):
#
setenv n  2
setenv m  1
setenv nm `echo $n $m | awk '{printf("%3.3d\n",$1*$2)}'`
echo "$n $m 0.75" | ../../../ALL/topo/src/partit
mv fort.21 depth_ATLb2.00_01.$nm
#
# --- 003 (1d):
#
setenv n  3
setenv m  1
setenv nm `echo $n $m | awk '{printf("%3.3d\n",$1*$2)}'`
echo "$n $m 0.75" | ../../../ALL/topo/src/partit
mv fort.21 depth_ATLb2.00_01.$nm
#
# --- 004:
#
setenv n  2
setenv m  2
setenv nm `echo $n $m | awk '{printf("%3.3d\n",$1*$2)}'`
echo "$n $m 0.75" | ../../../ALL/topo/src/partit
mv fort.21 depth_ATLb2.00_01.$nm
#
# --- 008:
#
setenv n  4
setenv m  2
setenv nm `echo $n $m | awk '{printf("%3.3d\n",$1*$2)}'`
echo "$n $m 0.75" | ../../../ALL/topo/src/partit
mv fort.21 depth_ATLb2.00_01.$nm
#
# --- 009:
#
setenv n  3
setenv m  3
setenv nm `echo $n $m | awk '{printf("%3.3d\n",$1*$2)}'`
echo "$n $m 0.75" | ../../../ALL/topo/src/partit
mv fort.21 depth_ATLb2.00_01.$nm
#
# --- 016:
#
setenv n  4
setenv m  4
setenv nm `echo $n $m | awk '{printf("%3.3d\n",$1*$2)}'`
echo "$n $m 0.75" | ../../../ALL/topo/src/partit
mv fort.21 depth_ATLb2.00_01.$nm
