#
set echo
#
# --- create list of sizes for dimensions.h.
#
head -1 depth_ATLb2.00_01.016 >! size.lis
#
foreach f ( *.*[0-9] )
  head -2 $f | tail -1 >> size.lis
end
#
#echo " " >> size.lis
##
#head -1 depth_ATLb2.00_01.016u >> size.lis
##
#foreach f ( *.*[0-9]u )
#  head -2 $f | tail -1 >> size.lis
#end
