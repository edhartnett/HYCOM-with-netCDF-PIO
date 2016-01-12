#!/bin/csh
#
#@ job_name         = 011w
#@ output           = $(job_name).log
#@ error            = $(job_name).log
#@ restart          = yes
#@ job_type         = serial
#@ node_usage       = shared
#@ wall_clock_limit = 12:00:00
#@ account_no       = NRLSS018
#@ notification     = never
#@ class            = share
#@ resources        = ConsumableCpus(1) ConsumableMemory(800mb)
#@ queue
#
set echo
set time = 1
set timestamp
C
C --- Create model interpolated winds for HYCOM. Use cubic cpline interpolation.
C --- Zero offset file.
C ---
C --- prebuild this script similar to a model script
# --- awk -f 284.awk y01=098 ab=a 284W.com > 284w098a.com
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
    if (-e ~wallcraf/bin/pget) then
      setenv pget ~wallcraf/bin/pget
      setenv pput ~wallcraf/bin/pput
    else
      setenv pget cp
      setenv pput cp
    endif
    breaksw
endsw
C
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C --- N is data-set name, e.g. ec10m-reanal
C --- W is permanent native wind directory
C
setenv E 020
setenv P hycom/ATLb2.00/expt_02.0/data
setenv D ~/$P
setenv N nogaps0.5c
#setenv N era40/6hourly
C
switch ($OS)
case 'SunOS':
    if (-e /net/hermes/scrb) then
#                  NRLSSC
      setenv S     /net/hermes/scrb/${user}/$P
      setenv W     /net/hermes/scrb/metzger/force/$N
    else
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv W     /u/home/metzger/force/$N
    endif
    breaksw
case 'Linux':
    if (-e /external/fast) then
#                  NRLSSC
      setenv S /external/fast/${user}/$P
      setenv W ~/force/$N
    else if (-e /work) then
#                  ERDC
      setenv S /work/${user}/$P
      setenv W ~/force/$N
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
      setenv W ~/force/$N
    endif
    breaksw
case 'OSF1':
    mkdir        ~/scratch
    chmod a+rx   ~/scratch
    setenv S     ~/scratch/$P
    setenv W     /u/home/metzger/force/$N
    breaksw
case 'IRIX64':
    mkdir        /workspace/${user}
    chmod a+rx   /workspace/${user}
    setenv S     /workspace/${user}/$P
    setenv W1    /msas031/metzger/force/$N1
    setenv W2    /msas031/metzger/force/$N2
    breaksw
case 'AIX':
    if (-e /gpfs/work) then
#                  ERDC MSRC
      mkdir        /gpfs/work/${user}
      chmod a+rx   /gpfs/work/${user}
      setenv S     /gpfs/work/${user}/$P
      setenv W     /u/home/metzger/force/$N
    else if (-e /scr) then
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv W     /u/home/metzger/force/$N
    else
#                  ARL MSRC, under GRD
      mkdir        /usr/var/tmp/${user}
      chmod a+rx   /usr/var/tmp/${user}
      setenv S     /usr/var/tmp/${user}/$P
      setenv W     /archive/navy/metzger/force/$N
    endif
    breaksw
case 'unicos':
    mkdir        /tmp/${user}
    chmod a+rx   /tmp/${user}
    mkdir        /tmp/${user}/ATLd0.08
    chgrp $ACCT  /tmp/${user}/ATLd0.08
    setenv S     /tmp/${user}/$P
    setenv W     /u/home/metzger/force/$N
    breaksw
endsw
C
mkdir -p $S/wind
cd       $S/wind
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
C --- For winds, only Y01 and A are used.
C
C --- One year spin-up run.
C
@ ymx =  1
C
setenv A "a"
setenv B "b"
setenv Y01 "099"
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
touch      fort.44 fort.44a fort.45 fort.45a
#/bin/rm -f fort.44 fort.44a fort.45 fort.45a
if (-z fort.44) then
  ${pget} $D/../../force/offset/tauewd_zero.b fort.44  &
endif
if (-z fort.44a) then
  ${pget} $D/../../force/offset/tauewd_zero.a fort.44a &
endif
if (-z fort.45) then
  ${pget} $D/../../force/offset/taunwd_zero.b fort.45  &
endif
if (-z fort.45a) then
  ${pget} $D/../../force/offset/taunwd_zero.a fort.45a &
endif
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
endif
C
@ i = 70
foreach y ( 2003 2004 2005 2006 2007 )
  @ i = $i + 1
  setenv N `echo $i | awk '{printf("%02d\n", $1)}'`
  touch  fort.${N}
  if (-z fort.${N}) then
    if     ($y >= 2003) then
      ${pget} ${W}/3hourly/nogaps0.5c-sea_${y}_03hr_strblk.D fort.${N} &
    else
      ${pget} ${W}/6hourly/nogaps1.0a-sea_${y}_06hr_strblk.D fort.${N} &
    endif
#   ${pget} ${W}/era40-sea_${y}_06hr_strblk.D fort.${N} &
  endif
end
C
C
C --- executable
C
/bin/cp /u/home/wallcraf/hycom/ALL/force/src/wi . &
wait
chmod ug+rx wi
ls -laFq
C
C --- NAMELIST input.
C
touch   fort.05in
/bin/rm fort.05in
cat <<E-o-D  > fort.05in
 &WWTITL
  CTITLE = '123456789012345678901234567890123456789012345678901234567890',
  CTITLE = 'NOGAPS 0.5c, 3-hrly, COARE v3.0, ocn-only, MKS',
 /
 &WWTIME
  SPDMIN =   0.0,  !minimum allowed wind speed
  WSCALE =   1.0,  !scale factor to mks
  WSTART = ${TS},
  TSTART = ${TS},
  TMAX   = ${TM},
 /
 &WWFLAG
  IGRID  = 2,  !0=p; 1=u&v; 2=p
  ISPEED = 0,  !0:none; 1:const; 2:kara; 3:coare
  INTERP = 1,  !0:bilinear; 1:cubic spline
  INTMSK = 0,  !0:no mask; 1:land/sea=0/1; 2:l/s=1/0;
  IFILL  = 3,  !0,1:tx&ty; 2,3:magnitude; 1,3:smooth; (intmsk>0 only)
  IOFILE = 0,  !0:single offset; 1:multiple offsets; 2:multi-off no .b check
  IWFILE = 4,  !1:ann/mon; 2:multi-file; 4:actual wind day
 /
E-o-D
switch ($OS)
case 'unicos':
case 'AIX':
C
C --- Fortran 90 NAMELIST delimiter
C
  /bin/rm -f fort.05i
  sed -e 's/&END/\//' -e 's/&end/\//' -e '/^#/d' < fort.05in > fort.05i
  breaksw
default:
C
C --- Fortran 77 NAMELIST delimiter
C
  /bin/rm -f fort.05i
  cp fort.05in fort.05i
  breaksw
endsw
C
C --- run the wind interpolation.
C
touch fort.10 fort.10a
/bin/rm -f fort.1[012] fort.1[012]a
C
setenv FOR010A fort.10a
setenv FOR011A fort.11a
setenv FOR012A fort.12a
C
setenv FOR044A fort.44a
setenv FOR045A fort.45a
C
switch ($OS)
case 'SunOS':
case 'Linux':
case 'OSF1':
case 'AIX':
    /bin/rm -f core
    touch core
    ./wi < fort.05i
    breaksw
case 'IRIX64':
    /bin/rm -f core
    touch core
    assign -V
    ./wi < fort.05i
    assign -V
    assign -R
    breaksw
case 'unicos':
    /bin/rm -f core
    touch core
    assign -V
    ./wi < fort.05i
    if (! -z core)  debug -s wi core
    assign -V
    assign -R
    breaksw
endsw
C
C --- Output.
C
/bin/mv fort.10  tauewd_${Y01}${A}.b
/bin/mv fort.10a tauewd_${Y01}${A}.a
/bin/mv fort.11  taunwd_${Y01}${A}.b
/bin/mv fort.11a taunwd_${Y01}${A}.a
if (-e fort.12) then
  /bin/mv fort.12  wndspd_${Y01}${A}.b
  /bin/mv fort.12a wndspd_${Y01}${A}.a
endif
C
if (-e ./SAVE) then
  ln tauewd_${Y01}${A}.a ./SAVE/tauewd_${Y01}${A}.a
  ln tauewd_${Y01}${A}.b ./SAVE/tauewd_${Y01}${A}.b
  ln taunwd_${Y01}${A}.a ./SAVE/taunwd_${Y01}${A}.a
  ln taunwd_${Y01}${A}.b ./SAVE/taunwd_${Y01}${A}.b
  if (-e wndspd_${Y01}${A}.a) then
    ln wndspd_${Y01}${A}.a ./SAVE/wndspd_${Y01}${A}.a
    ln wndspd_${Y01}${A}.b ./SAVE/wndspd_${Y01}${A}.b
  endif
endif
C
C  --- END OF JUST IN TIME WIND GENERATION SCRIPT.
C
#
/usr/lpp/LoadL/full/bin/llq -w $LOADL_STEP_ID
