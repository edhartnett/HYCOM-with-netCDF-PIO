#
set echo
#
# --- build new topo scripts from old.
# --- foreach should be customized for needed scripts.
# --- some scripts will need additional manual editing.
#
# RO = old region
# RN = new region
#
setenv RO ATLa2.00
setenv RN ATLb2.00
#
cd ../../${RO}/topo
#
sed -e "s/${RO}/${RN}/g" regional.grid.com >! ../../${RN}/topo/regional.grid.com
#
foreach t ( _01.com _02.com _03.com )
  sed -e "s/${RO}/${RN}/g" depth_${RO}$t >! ../../${RN}/topo/depth_${RN}$t
end
#
#foreach t ( .com _modify.com )
#  sed -e "s/${RO}/${RN}/g" landsea_${RO}$t >! ../../${RN}/topo/landsea_${RN}$t
#end
