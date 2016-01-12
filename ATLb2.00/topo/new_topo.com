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
#
cd partit
mkdir ../../../${RN}/topo/partit
#
foreach t ( _2d.com )
  sed -e "s/${RO}/${RN}/g" -e "s/_01/_01/g" depth_${RO}_01$t >! ../../../${RN}/topo/partit/depth_${RN}_01$t
end
foreach f ( ppm.com size.com )
  sed -e "s/${RO}/${RN}/g" -e "s/_01/_01/g" $f >! ../../../${RN}/topo/partit/$f
end
cp xbathy.pal ../../../${RN}/topo/partit/xbathy.pal
