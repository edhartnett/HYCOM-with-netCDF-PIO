#!/bin/csh
set echo
#
# --- build new expt files from old.
# --- some files will need additional manual editing.
#
# DO = old experiment directory name
#  O = old experiment number
# DN = new experiment directory name
#  N = new experiment number
#  R = region name.
#
setenv DO expt_01.3
setenv  O 013
setenv DN expt_01.5
setenv  N 015
setenv  R ATLb2.00
#
setenv RO $R
setenv  D ../../${RO}/${DO}

foreach t ( .com W.com F.com P.com y001.limits cod.com lsf.com nqs.com rll.com grd.com )
  sed -e "s/setenv E .*${O}.*/setenv E ${N}/g" -e "s/${DO}/${DN}/g"  -e "s/${RO}/${R}/g" ${D}/${O}${t} >! ${N}${t}
end
#
cp ${D}/${O}.awk ${N}.awk
#
cp ${D}/blkdat.input ${D}/ports.input* .
#
cp ${D}/dummy*.com .
#
# --- experiment run sequence:
#
mlist 001 020 1
cp LIST LIST++
#
# --- create data directory (local and on archive):
#
mkdir data
#krsh vincent.navo.hpc.mil mkdir -p /u/home/$user/hycom/${R}/${DN}/data
