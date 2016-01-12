#!/bin/csh
#
set echo
set time = 1
set timestamp
C
C --- Create model interpolated winds for HYCOM.
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
case 'IRIX64':
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
default:
    setenv pget cp
    setenv pput cp
endsw
C
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C --- X is data-set executable abbreviation, e.g. 1125_ec
C --- N is data-set name, e.g. ec10m-reanal
C --- W is permanent native pcip directory
C
setenv E 015
setenv P hycom/ATLb2.00/expt_01.5/data
setenv D ~/$P
setenv X 1125_ec
setenv N ec10m-reanal
C
switch ($OS)
case 'SunOS':
    if (-e /net/hermes/scrb) then
#                  NRLSSC
      setenv S     /net/hermes/scrb/${user}/$P
      setenv W     /net/hermes/scrb/metzger/wind_ieee/$N
    else
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv W     /u/b/metzger/wind_ieee/$N
    endif
    breaksw
case 'Linux':
    if (-e /external/fast) then
#                  NRLSSC
      setenv S /external/fast/${user}/$P
      setenv W ~/wind_ieee/$N
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
      setenv W ~/wind_ieee/$N
    endif
    breaksw
case 'OSF1':
    mkdir        ~/scratch
    chmod a+rx   ~/scratch
    setenv S     ~/scratch/$P
    setenv W     /u/b/metzger/wind_ieee/$N
    breaksw
case 'IRIX64':
    mkdir        /scr/${user}
    chmod a+rx   /scr/${user}
    setenv S     /scr/${user}/$P
    setenv W     /u/b/metzger/wind_ieee/$N
    breaksw
case 'AIX':
    if (-e /gpfs/work) then
#                  ERDC MSRC
      mkdir        /gpfs/work/${user}
      chmod a+rx   /gpfs/work/${user}
      setenv S     /gpfs/work/${user}/$P
      setenv W     /u/b/metzger/wind_ieee/$N
    else if (-e /scr) then
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv W     /u/b/metzger/wind_ieee/$N
    else
#                  ARL MSRC, under GRD
      mkdir        /usr/var/tmp/${user}
      chmod a+rx   /usr/var/tmp/${user}
      setenv S     /usr/var/tmp/${user}/$P
      setenv W     /u/b/metzger/wind_ieee/$N
    endif
    breaksw
case 'unicos':
case 'unicosmk':
    mkdir        /tmp/${user}
    chmod a+rx   /tmp/${user}
    mkdir        /tmp/${user}/ATLd0.08
    chgrp $ACCT  /tmp/${user}/ATLd0.08
    setenv S     /tmp/${user}/$P
    setenv W     /u/b/metzger/wind_ieee/$N
    breaksw
endsw
C
mkdir -p $S/wind
cd       $S/wind
C
C --- For whole year runs.
C ---   Y01 initial model year of this run.
C ---   YXX is the last model year of this run, and the first of the next run.
C ---   A and B are identical, typically blank.
C --- For part year runs.
C ---   A is this part of the year, B is next part of the year.
C ---   Y01 is the start model year of this run.
C ---   YXX is the end   model year of this run, usually Y01.
C --- For a few hour/day run
C ---   A   is the start day and hour, of form "dDDDhHH".
C ---   B   is the end   day and hour, of form "dXXXhYY".
C ---   Y01 is the start model year of this run.
C ---   YXX is the end   model year of this run, usually Y01.
C --- Note that these variables are set by the .awk generating script.
C
C --- For winds, only Y01 and A are used.
C
setenv A "a"
setenv B "b"
setenv Y01 "099"
setenv YXX "099"
C
echo "Y01 =" $Y01 "YXX = " $YXX  "A =" ${A} "B =" ${B}
C
C --- time limits.
C
if (-e ${D}/../${E}y${Y01}${A}.limits) then
  setenv TS `cat ${D}/../${E}y${Y01}${A}.limits | awk '{print $1}'`
  setenv TM `cat ${D}/../${E}y${Y01}${A}.limits | awk '{print $2}'`
else
# use "LIMITI" when starting a run after day zero.
# use "LIMITS9" (say) for a 9-day run.
  setenv TS `echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} | awk '{print $1}'`
  setenv TM `echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} | awk '{print $2}'`
endif
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
touch      fort.71 fort.72 fort.73 fort.74 fort.75 fort.76
if (-z fort.71) then
  ${pget} $W/ds115.3_ec10m-reanal_7982_06hr_ws.d fort.71 &
endif
if (-z fort.72) then
  ${pget} $W/ds115.3_ec10m-reanal_8386_06hr_ws.d fort.72 &
endif
if (-z fort.73) then
  ${pget} $W/ds115.3_ec10m-reanal_8790_06hr_ws.d fort.73 &
endif
if (-z fort.74) then
  ${pget} $W/ds115.3_ec10m-reanal_9193_06hr_ws.d fort.74 &
endif
if (-z fort.75) then
  ${pget} $W/../ec10m-oper/ds111.1_ec10m-oper_9497_06hr_ws.d fort.75 &
endif
if (-z fort.76) then
  ${pget} $W/../ec10m-oper/ds111.1_ec10m-oper_9801_06hr_ws.d fort.76 &
endif
C
C --- executable
C
####${pget} ${D}/../../../ALL/force/src/wi_${X} . &
/bin/cp ${D}/../../../ALL/force/src/wi_${X} . &
wait
chmod ug+rx wi_${X}
ls -laFq
C
C --- NAMELIST input.
C
touch   fort.05in
/bin/rm fort.05in
cat <<E-o-D  > fort.05in
 &WWTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = 'ECMWF-Reanalysis, 10m DS111.1, 6-hr, MKS',
 &END
 &WWTIME
  SPDMIN = 0.0,
  WSTART = ${TS},
  TSTART = ${TS},
  TMAX   = ${TM},
 &END
 &WWFLAG
  INTERP = 1
  IWFILE = 4,
  IWTILE = 1,
 &END
E-o-D
switch ($OS)
case 'unicos':
case 'unicos-lanl':
case 'unicosmk':
case 'unicos-t3d':
case 'sn6705':
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
    ./wi_${X} < fort.05i
    breaksw
case 'IRIX64':
    /bin/rm -f core
    touch core
    assign -V
    ./wi_${X} < fort.05i
    assign -V
    assign -R
    breaksw
case 'unicosmk':
    /bin/rm -f core
    touch core
    assign -V
    ./wi_${X} < fort.05i
    if (! -z core)  debugview wi_${X} core
    assign -V
    assign -R
    breaksw
case 'unicos':
    /bin/rm -f core
    touch core
    assign -V
    ./wi_${X} < fort.05i
    if (! -z core)  debug -s wi_${X} core
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
/bin/mv fort.12  wndspd_${Y01}${A}.b
/bin/mv fort.12a wndspd_${Y01}${A}.a
C
if (-e ./SAVE) then
  ln tauewd_${Y01}${A}.a ./SAVE/tauewd_${Y01}${A}.a
  ln tauewd_${Y01}${A}.b ./SAVE/tauewd_${Y01}${A}.b
  ln taunwd_${Y01}${A}.a ./SAVE/taunwd_${Y01}${A}.a
  ln taunwd_${Y01}${A}.b ./SAVE/taunwd_${Y01}${A}.b
  ln wndspd_${Y01}${A}.a ./SAVE/wndspd_${Y01}${A}.a
  ln wndspd_${Y01}${A}.b ./SAVE/wndspd_${Y01}${A}.b
endif
C
C  --- END OF JUST IN TIME WIND GENERATION SCRIPT.
C
