#! /bin/csh
#BSUB -n 1
#BSUB -W 24:00
set echo
#
# --- form the mean and mean-square of a sequence of archive files
#
# --- R is region name.
# --- T is topography number.
# --- K is number of layers (or 1 for a surface mean).
# --- E is expt number.
# --- P is primary path.
# --- D is permanent directory.
# --- S is scratch   directory, can be the permanent directory.
#
setenv R ATLb2.00
setenv T 01
setenv K 22
setenv E 010
setenv P hycom/${R}/expt_01.0/data
setenv D ~/$P
setenv S $D
#
# --- pget, pput "copy" files between scratch and permanent storage.
#
setenv pget "/bin/ln -s"
setenv pput "/bin/cp"
#
mkdir -p  ${S}
cd        ${S}
#
cp ${D}/../${E}.awk .
#
touch   regional.depth.a regional.depth.b
/bin/rm regional.depth.a regional.depth.b
ln -s ${D}/../../topo/depth_${R}_${T}.a regional.depth.a
ln -s ${D}/../../topo/depth_${R}_${T}.b regional.depth.b
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
# --- individual model years, 30-day 3-d archives.
#
setenv YRFLAG 0
setenv ARCHFQ 30.0
#
foreach Y ( 0020 )
#foreach mon ( a b c d e f g h i j k l )
setenv mon ""
  echo "month " ${mon}
#
  setenv DF `echo LIMITS | awk -f ${E}.awk y01=${Y} ab=${mon} | awk '{print $1 - 0.01}'`
  setenv DL `echo LIMITS | awk -f ${E}.awk y01=${Y} ab=${mon} | awk '{print $2 - 0.01}'`
#
  echo $YRFLAG $ARCHFQ $DF $DL | ${D}/../../../ALL/bin/hycom_archv_dates | tail +3 >! ${E}_mean.tmp
  foreach d ( `cat ${E}_mean.tmp` )
    touch archv.${d}.a
    if (-z archv.${d}.a) then
      /bin/rm archv.${d}.a
      ${pget} ${D}/archv.${d}.a archv.${d}.a &
    endif
    touch archv.${d}.b
    if (-z archv.${d}.b) then
      /bin/rm archv.${d}.b
      ${pget} ${D}/archv.${d}.b archv.${d}.b &
    endif
  end
  wait
#
# --- mean.
#
  setenv NA `cat ${E}_mean.tmp | wc -l`
  echo "  $K	'kk    ' = number of layers involved"               >! ${E}_mean.IN
  echo "  0	'meansq' = form meansq, rather than mean (0=F,1=T)" >> ${E}_mean.IN
  echo "  $NA	'narchs' = number of archives to read (==0 to end input)" >> ${E}_mean.IN
  cat ${E}_mean.tmp | awk '{printf("archv.%s.a\n",$1)}' >> ${E}_mean.IN
  echo "  0  	'narchs' = number of archives to read (==0 to end input)" >> ${E}_mean.IN
  echo   "archMN.${Y}${mon}" >> ${E}_mean.IN
  touch   archMN.${Y}${mon}.a archMN.${Y}${mon}.b
  /bin/rm archMN.${Y}${mon}.a archMN.${Y}${mon}.b
  cat ${E}_mean.IN
  ${D}/../../../ALL/meanstd/src/hycom_mean < ${E}_mean.IN
  cat ${E}_mean.IN
  if ($D != $S) then
    ${pput} archMN.${Y}.a ${D}/archMN.${Y}.a
    ${pput} archMN.${Y}.b ${D}/archMN.${Y}.b
  endif
#
# --- mnsq.
#
  sed -e "s/0	'meansq'/1	'meansq'/g" -e "s/archMN/archSQ/g" ${E}_mean.IN >! ${E}_mnsq.IN
  touch   archSQ.${Y}${mon}.a archSQ.${Y}${mon}.b
  /bin/rm archSQ.${Y}${mon}.a archSQ.${Y}${mon}.b
  cat ${E}_mnsq.IN
  ${D}/../../../ALL/meanstd/src/hycom_mean < ${E}_mnsq.IN
  cat ${E}_mnsq.IN
  if ($D != $S) then
    ${pput} archSQ.${Y}.a ${D}/archSQ.${Y}.a
    ${pput} archSQ.${Y}.b ${D}/archSQ.${Y}.b
  endif
#
# /bin/rm ${E}_mean.IN ${E}_mnsq.IN ${E}_mean.tmp
#end
end
