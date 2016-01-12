#!/bin/csh
#
set echo
set time = 1
set timestamp
C
C --- Create model interpolated surtmp for HYCOM.
C ---
C --- prebuild this script similar to a model script
# --- awk -f 284.awk y01=098 ab=a 284T.com > 284t098a.com
C
C --- Preamble, script keys on O/S name.
C
setenv OS `uname`
switch ($OS)
case 'SunOS':
C   assumes /usr/5bin is before /bin and /usr/bin in PATH.
    breaksw
case 'Linux':
    breaksw
case 'OSF1':
    breaksw
case 'IRIX64':
    breaksw
case 'AIX':
    breaksw
case 'unicosmk':
    setenv ACCT `groups | awk '{print $1}'`
    breaksw
case 'unicos':
    setenv ACCT `newacct -l | awk '{print $4}'`
    breaksw
default:
    echo 'Unknown Operating System: ' $OS
    exit (1)
endsw
C
C --- pget, pput "copy" files between scratch and permanent storage.
C --- Can both be cp if the permanent filesystem is mounted locally.
C
switch ($OS)
case 'SunOS':
case 'Linux':
case 'OSF1':
case 'AIX':
case 'unicos':
case 'unicosmk':
    if (-e ~wallcraf/bin/pget) then
      setenv pget ~wallcraf/bin/pget
      setenv pput ~wallcraf/bin/pput
    else
      setenv pget cp
      setenv pput cp
    endif
    breaksw
case 'IRIX64':
    setenv pget cp
    setenv pput cp
    breaksw
default:
    setenv pget cp
    setenv pput cp
endsw
C
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C --- N is data-set name, e.g. ecmwf-reanal_ds111.6
C --- W is permanent native surtmp directory
C
setenv E 016
setenv P hycom/ATLb2.00/expt_01.6/data
setenv D ~/$P
setenv N ecmwf-oper_ds111.1
C
switch ($OS)
case 'SunOS':
    if (-e /net/hermes/scrb) then
#                  NRLSSC
      setenv S     /net/hermes/scrb/${user}/$P
      setenv W     /net/hermes/scrb/metzger/temp_ieee/$N
    else
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv W     /u/home/metzger/temp_ieee/$N
    endif
    breaksw
case 'Linux':
    if (-e /external/fast) then
#                  NRLSSC
      setenv S /external/fast/${user}/$P
      setenv W ~/temp_ieee/$N
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
      setenv W ~/temp_ieee/$N
    endif
    breaksw
case 'OSF1':
#                 ERDC MSRC
    mkdir        /work/${user}
    chmod a+rx   /work/${user}
    setenv S     /work/${user}/$P
    setenv W     /u/home/metzger/temp_ieee/$N
    breaksw
case 'IRIX64':
    mkdir        /workspace/${user}
    chmod a+rx   /workspace/${user}
    setenv S     /workspace/${user}/$P
    setenv W     /msas031/metzger/temp_ieee/$N
    breaksw
case 'AIX':
    if (-e /gpfs/work) then
#                  ERDC MSRC
      mkdir        /gpfs/work/${user}
      chmod a+rx   /gpfs/work/${user}
      setenv S     /gpfs/work/${user}/$P
      setenv W     /u/home/metzger/temp_ieee/$N
    else if (-e /scr) then
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv W     /u/home/metzger/temp_ieee/$N
    else
#                  ARL MSRC, under GRD
      mkdir        /usr/var/tmp/${user}
      chmod a+rx   /usr/var/tmp/${user}
      setenv S     /usr/var/tmp/${user}/$P
      setenv W     /u/home/metzger/temp_ieee/$N
    endif
    breaksw
case 'unicos':
case 'unicosmk':
    mkdir        /tmp/${user}
    chmod a+rx   /tmp/${user}
    mkdir        /tmp/${user}/ATLb2.00
    chgrp $ACCT  /tmp/${user}/ATLb2.00
    setenv S     /tmp/${user}/$P
    setenv W     /u/home/metzger/temp_ieee/$N
    breaksw
endsw
C
mkdir -p $S/ssta
cd       $S/ssta
C
C --- For whole year runs.
C ---   ymx number of years per model run.
C ---   Y01 initial model year of this run.
C ---   YXX is the last model year of this run, and the first of the next run.
C ---   A and B are identical, typically blank.
C --- For part year runs.
C ---   A is this part of the year, B is next part of the year.
C ---   Y01 initial model year of this run.
C ---   YXX year at end of this part year run.
C ---   ymx is 1.
C --- Note that these variables and the .awk generating script must
C ---  be consistant to get a good run script.
C
C --- For surtmp, only Y01 and A are used.
C
C --- One year spin-up run.
C
@ ymx =  1
C
setenv A "a"
setenv B "b"
setenv Y01 "096"
C
switch ("${B}")
case "${A}":
    setenv YXX `echo $Y01 $ymx | awk '{printf("%03d", $1+$2)}'`
    breaksw
case "a":
    setenv YXX `echo $Y01 | awk '{printf("%03d", $1+1)}'`
    breaksw
default:
    setenv YXX $Y01
endsw
C
echo "Y01 =" $Y01 "YXX = " $YXX  "A =" ${A} "B =" ${B}
C
C --- time limits.
C --- use "LIMITI" when starting a run after day zero.
C
setenv TS `echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} | awk '{print $1}'`
setenv TM `echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} | awk '{print $2}'`
C
echo "TS =" $TS "TM =" $TM
C
C --- input files from file server.
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
endif
C
touch      fort.70
if (-z fort.70) then
  ${pget} ${W}/../../mask_ieee/ecmwf-op/ecmwf-oper_landsea_T319.D fort.70 &
endif
C
touch  fort.71 fort.72 
if (-z fort.71) then
  ${pget} ${W}/ds111.1_9801_06hr_soiltemp.d fort.71 &
endif
if (-z fort.72) then
  ${pget} ${W}/ds111.1_0205_06hr_soiltemp.d fort.72 &
endif
C
C --- executable
C
/bin/cp /u/home/wallcraf/hycom/ALL/force/src/kp . &
wait
chmod ug+rx kp
ls -laFq
C
C --- NAMELIST input.
C
cat <<E-o-D  >! fort.05i
 &AFTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = 'ECMWF soiltmp, masked, 6hrly, degC',
  CNAME  = 'surtmp',
 /
 &AFTIME
  PARMIN = -999.0,  ! disable parmin
  PARMAX =  999.0,  ! disable parmax
  PAROFF = -273.15, ! Kelvins to degC
  TMPICE =   -1.7,  !ECMWF sea ice marker
  FSTART = ${TS},
  TSTART = ${TS},
  TMAX   = ${TM},
 /
 &AFFLAG
  IFFILE =   5,
  INTERP =   0,  !0:piecewise-linear; 1:cubic-spline;
  INTMSK =   2,  !0:no mask; 1:land/sea=0/1; 2:land/sea=1/0;
 /
E-o-D
C
C --- run the surtmp interpolation.
C
touch fort.10 fort.10a
/bin/rm -f fort.10 fort.10a
C
setenv FOR010A fort.10a
C
switch ($OS)
case 'SunOS':
case 'Linux':
case 'OSF1':
case 'AIX':
    /bin/rm -f core
    touch core
    ./kp < fort.05i
    breaksw
case 'IRIX64':
    /bin/rm -f core
    touch core
    assign -V
    ./kp < fort.05i
    assign -V
    assign -R
    breaksw
case 'unicosmk':
    /bin/rm -f core
    touch core
    assign -V
    ./kp < fort.05i
    if (! -z core)  debugview kp core
    assign -V
    assign -R
    breaksw
case 'unicos':
    /bin/rm -f core
    touch core
    assign -V
    ./kp < fort.05i
    if (! -z core)  debug -s kp core
    assign -V
    assign -R
    breaksw
endsw
C
C --- Output.
C
/bin/mv fort.10  surtmp_${Y01}${A}.b
/bin/mv fort.10a surtmp_${Y01}${A}.a
C
if (-e ./SAVE) then
  ln surtmp_${Y01}${A}.a ./SAVE/surtmp_${Y01}${A}.a
  ln surtmp_${Y01}${A}.b ./SAVE/surtmp_${Y01}${A}.b
endif
C
C  --- END OF JUST IN TIME SURTMP GENERATION SCRIPT.
C
