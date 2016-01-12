#
set echo
#
foreach f ( READ* Makefile *.h *.f *.F )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -ibw $f ../src_2.0.01_16_omp
end
