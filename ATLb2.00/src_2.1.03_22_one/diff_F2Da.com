#
set echo
#
foreach f ( READ* Makefile *.h *.f *.F )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -ibw $f ../../F2Da0.03/src_2.0.03_50_one
end
