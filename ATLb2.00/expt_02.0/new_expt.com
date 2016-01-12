#!/bin/csh
set echo
#
# --- build new expt files from old.
# --- some files will need additional manual editing.
#
#  O = old experiment number
#  N = new experiment number
#  R = region name.
#
setenv  O 020
setenv  N 021
setenv  R ATLb2.00
#
setenv DO `echo ${O} | awk '{printf("expt_%04.1f", $1*0.1)}'`
setenv DN `echo ${N} | awk '{printf("expt_%04.1f", $1*0.1)}'`
#
setenv RO $R
setenv  D ../../${RO}/${DO}
#
foreach t ( .com A.com F.com O.com P.com Q.com T.com W.com S.com y001.limits cod.com lsf.com nqs.com rll.com grd.com pbs.com )
  sed -e "s/setenv E .*${O}.*/setenv E ${N}/g" -e "s/${DO}/${DN}/g"  -e "s/${RO}/${R}/g" ${D}/${O}${t} >! ${N}${t}
end
#
cp ${D}/${O}.awk ${N}.awk
#
cp ${D}/blkdat.input* ${D}/ports.input* ${D}/tracer.input* .
#
cp ${D}/dummy*.com .
#
# --- experiment run sequence:
#
mlist 001 020 1 
#cp  ${D}/LIST++ LIST
cp LIST LIST++
#
# --- create data directory (local and on archive):
#
mkdir data
#krsh newton.navo.hpc.mil mkdir -p /u/home/$user/hycom/${R}/${DN}/data
