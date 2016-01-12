#!/bin/csh
set echo
#
# --- build new atmospheric forcing scripts from old.
# --- foreach should be customized for needed scripts.
# --- some scripts will need additional manual editing.
#
# RO = old region
# RN = new region
#
setenv RO ATLb2.00
setenv RN GLBx2.00
#
mkdir offset
mkdir coads
mkdir plot
#
cd ../../${RO}/force
#
foreach f ( offset/tauXwd_zero.com coads/precip_zero.com coads/coads_mon_wind.com coads/coads_mon_flux.com )
  sed -e "s/${RO}/${RN}/g" $f >! ../../${RN}/force/$f
end
#
# --- change idm and jdm sizes for these particular cases.
#
foreach f ( plot/link.com plot/alias.src plot/*.IN )
  sed -e "s/${RO}/${RN}/g" -e "s/ 57	'idm   '/180	'idm   '/g" -e "s/ 52	'jdm   '/176	'jdm   '/g" -e "s/10	'lalolb'/20	'lalolb'/g" -e "s/ 10	'lalogr'/-10	'lalogr'/g" -e "s/4	'loclab'/3	'loclab'/g" -e "s/11	'locbar'/11	'locbar'/g" $f >! ../../${RN}/force/$f
end
