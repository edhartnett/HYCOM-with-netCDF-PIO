#! /bin/csh
#BSUB -n 1
#BSUB -W 24:00
set echo
#
# --- form the mean and mean-square of a sequence of archive files
#
# --- R is region name.
# --- K is number of layers (or 1 for a surface mean).
# --- P is primary path.
# --- D is permanent directory.
# --- S is scratch   directory, can be the permanent directory.
# --- Y is mean extent.
#
setenv R ATLb2.00
setenv K 22
setenv P hycom/${R}/expt_01.0/data
setenv D ~/$P
setenv S $D
setenv Y 0020
#
# --- pget, pput "copy" files between scratch and permanent storage.
#
setenv pget "/bin/ln -s"
setenv pput "/bin/cp"
#
mkdir -p  ${S}
cd        ${S}
#
touch regional.grid.a regional.grid.b
if (-z regional.grid.a) then
  /bin/rm regional.grid.a
  /bin/ln -s ${D}/../../topo/regional.grid.a regional.grid.a
endif
if (-z regional.grid.b) then
  /bin/rm regional.grid.b
  /bin/ln -s ${D}/../../topo/regional.grid.b regional.grid.b
endif
#
touch archMN.${Y}.a archMN.${Y}.b archSQ.${Y}.a archSQ.${Y}.b
if (-z archMN.${Y}.a) then
  /bin/rm archMN.${Y}.a
  ${pget} ${D}/archMN.${Y}.a archMN.${Y}.a &
endif
if (-z archMN.${Y}.b) then
  /bin/rm archMN.${Y}.b
  ${pget} ${D}/archMN.${Y}.b archMN.${Y}.b &
endif
if (-z archSQ.${Y}.a) then
  /bin/rm archSQ.${Y}.a
  ${pget} ${D}/archSQ.${Y}.a archSQ.${Y}.a &
endif
if (-z archSQ.${Y}.b) then
  /bin/rm archSQ.${Y}.b
  ${pget} ${D}/archSQ.${Y}.b archSQ.${Y}.b &
endif
wait
#
touch   ${D}/archSD.${Y}.a ${D}/archSD.${Y}.b
/bin/rm ${D}/archSD.${Y}.a ${D}/archSD.${Y}.b
${D}/../../../ALL/meanstd/src/hycom_std << E-o-D
 $K	'kk    ' = number of layers involved
archMN.${Y}.b
archSQ.${Y}.b
archSD.${Y}
E-o-D
#
if ($D != $S) then
  ${pput} archSD.${Y}.a ${D}/archSD.${Y}.a
  ${pput} archSD.${Y}.b ${D}/archSD.${Y}.b
endif
