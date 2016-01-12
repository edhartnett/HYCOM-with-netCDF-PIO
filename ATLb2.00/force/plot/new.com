#!/bin/csh
#
set echo
#
# --- build new plot scripts from old.
# --- foreach should be customized for needed scripts.
# --- some scripts will need additional manual editing.
#
foreach f ( *.IN )
  /bin/mv ${f} ${f}+
  sed    -e "s/ATLb2.00/PORd0.32/g" \
      -e "s/^.*'idm   '/ 69	'idm   '/" \
      -e "s/^.*'jdm   '/ 69	'jdm   '/" \
      -e "s/^.*'lalolb'/  2	'lalolb'/" \
      -e "s/^.*'lalogr'/ -2	'lalogr'/" \
      -e "s/^.*'loclab'/  2	'loclab'/" \
      -e "s/^.*'locbar'/ 14	'locbar'/" \
      ${f}+ > ${f}
  /bin/rm ${f}+
end
