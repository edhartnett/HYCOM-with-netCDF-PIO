#
set echo
#
# --- form the mean from a HYCOM transport section file
#
setenv D  ../expt_01.0/data
setenv E  010
setenv YF 0020
setenv YL 0020
#
# --- standard "all layer" means
#
touch   ${D}/trans.${YF}-${YL}.mn
/bin/rm ${D}/trans.${YF}-${YL}.mn
./src/meantspt <<E-o-D
${D}/trans.${YF}-${YL}.d
${D}/trans.${YF}-${YL}.mn
  22	  'kout  ' = number of layer combinations to output (<=kdm)
E-o-D
#
# --- combine layers into 5 ranges
#
touch   ${D}/trans.${YF}-${YL}.mn05
/bin/rm ${D}/trans.${YF}-${YL}.mn05
./src/meantspt <<E-o-D
${D}/trans.${YF}-${YL}.d
${D}/trans.${YF}-${YL}.mn05
   5	  'kout  ' = number of layer combinations to output (<=kdm)
   8	  'laybot' = last layer in next combination of layers
  11	  'laybot' = last layer in next combination of layers
  15	  'laybot' = last layer in next combination of layers
  18	  'laybot' = last layer in next combination of layers
  22	  'laybot' = last layer in next combination of layers
E-o-D
