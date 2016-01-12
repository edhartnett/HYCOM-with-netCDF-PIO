#!/bin/csh
#
#set echo
#
# --- create all plots.
#
foreach n ( coads )
  foreach t ( airtmp vapmix )
    env NCARG_GKS_PS=${n}_${t}.ps ./src/fieldproc < ${n}_${t}.IN >&! ${n}_${t}.log &
  end
  wait
  foreach t ( radflx shwflx )
    env NCARG_GKS_PS=${n}_${t}.ps ./src/fieldproc < ${n}_${t}.IN >&! ${n}_${t}.log &
  end
  wait
  foreach t ( precip wndspd )
    env NCARG_GKS_PS=${n}_${t}.ps ./src/fieldproc < ${n}_${t}.IN >&! ${n}_${t}.log &
  end
  wait
end
