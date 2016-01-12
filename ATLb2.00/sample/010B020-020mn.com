#
set echo
#
# --- form the mean from a HYCOM MB transport section file
# --- these only contain the barotropic transport
#
setenv D  ../expt_01.0/data
setenv E  010
setenv Y 0020
#
# --- standard "all layer" means
#
touch   ${D}/transMB.${Y}.mn
/bin/rm ${D}/transMB.${Y}.mn
./src/meantspt <<E-o-D
${D}/transMB.${Y}.d
${D}/transMB.${Y}.mn
  22	  'kout  ' = number of layer combinations to output (<=kdm)
E-o-D
#
# --- smallest possible mean file
#
touch   ${D}/transMB.${Y}.mn03
/bin/rm ${D}/transMB.${Y}.mn03
./src/meantspt <<E-o-D
${D}/transMB.${Y}.d
${D}/transMB.${Y}.mn03
   3	  'kout  ' = number of layer combinations to output (<=kdm)
   1	  'laybot' = last layer in next combination of layers
   2	  'laybot' = last layer in next combination of layers
  22	  'laybot' = last layer in next combination of layers
E-o-D
