#
set echo
#
# --- form the mean from a HYCOM MN transport section file
#
setenv D  ../expt_01.0/data
setenv E  010
setenv Y 0020
#
# --- standard "all layer" means
#
touch   ${D}/transMN.${Y}.mn
/bin/rm ${D}/transMN.${Y}.mn
./src/meantspt <<E-o-D
${D}/transMN.${Y}.d
${D}/transMN.${Y}.mn
  22	  'kout  ' = number of layer combinations to output (<=kdm)
E-o-D
#
# --- combine layers into 5 ranges
#
touch   ${D}/transMN.${Y}.mn05
/bin/rm ${D}/transMN.${Y}.mn05
./src/meantspt <<E-o-D
${D}/transMN.${Y}.d
${D}/transMN.${Y}.mn05
   5	  'kout  ' = number of layer combinations to output (<=kdm)
   8	  'laybot' = last layer in next combination of layers
  11	  'laybot' = last layer in next combination of layers
  15	  'laybot' = last layer in next combination of layers
  18	  'laybot' = last layer in next combination of layers
  22	  'laybot' = last layer in next combination of layers
E-o-D
