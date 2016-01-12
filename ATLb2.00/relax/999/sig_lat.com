#
#set echo
#
# --- extract values needed for poflat
# --- uniform latitude spacing
#
echo "********************"
@ l = -30
while ($l < 66)
  foreach f ( sig_2[1-7].[05] )
    awk -f sig_lat.awk lat=$l $f
  end
  echo $l | awk '{printf(" 28.0 %6.1f 6000.0\n",$1)}'
  echo "********************"
  @ l = $l + 5
end
