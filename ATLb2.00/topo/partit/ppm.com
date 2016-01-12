#
set echo
#
# --- generate full size PPM image file.
#
cd ~/hycom/ATLb2.00/topo/partit
#
setenv FOR021  fort.21
setenv FOR051  ../depth_ATLb2.00_01.b
setenv FOR051A ../depth_ATLb2.00_01.a
foreach f ( *[0-9]? )
  cp $f fort.21
# cat fort.21
  ../../../ALL/topo/src/topo_ppm
  mv fort.31 ${f}.ppm
  /bin/rm -f fort.21
end
