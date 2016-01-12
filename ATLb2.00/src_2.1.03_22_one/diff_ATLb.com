#
set echo
#
foreach f ( READ* Makefile *.h *.f *.F )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -ibw $f ../../ATLb2.00/src_2.1.02_22_one
end
