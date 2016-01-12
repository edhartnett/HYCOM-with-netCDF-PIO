#
set echo
#
# --- do all the 010 plots.
#
foreach t ( hor cs1 cs2 )
  foreach p ( 010_jan 010_jul )
    env NCARG_GKS_PS=${p}_${t}.ps ./src/hycomproc < ${p}_${t}.IN >&! ${p}_${t}.log &
  end
  wait
end
