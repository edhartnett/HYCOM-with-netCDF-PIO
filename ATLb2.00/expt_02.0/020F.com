#!/bin/csh
#
#@ job_name         = 011f
#@ output           = $(job_name).log
#@ error            = $(job_name).log
#@ restart          = yes
#@ job_type         = serial
#@ node_usage       = shared
#@ wall_clock_limit = 12:00:00
#@ account_no       = NRLSS018
#@ notification     = never
#@ class            = share
#@ resources        = ConsumableCpus(1) ConsumableMemory(750mb)
#@ queue
#
set echo
set time = 1
set timestamp
C
C --- Create model interpolated fluxes for HYCOM. Use linear interpolation
C --- on fluxes because they change on such short space scales.
C
C --- Also CICE interpolated fluxes, if needed.  Note that these include
C --- idm and jdm, and so need customizing for the model domain.
C ---
C --- offset annual mean airtmp, radflx, shwflx.
C ---
C --- prebuild this script similar to a model script
# --- awk -f 284.awk y01=098 ab=a 284F.com > 284f098a.com
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
default:
    setenv pget cp
    setenv pput cp
endsw
C
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C --- N is data-set name, e.g. ecmwf-reanal
C --- W is permanent native flux directory
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
mkdir -p $S/flux
cd       $S/flux
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
C --- For fluxes, only Y01 and A are used.
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
C
if (-e ${D}/../${E}y${Y01}${A}.limits) then
  setenv TS `sed -e "s/-/ /g" ${D}/../${E}y${Y01}${A}.limits | awk '{print $1}'`
  setenv TM `cat              ${D}/../${E}y${Y01}${A}.limits | awk '{print $2}'`
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
      ${pget} ${W}/3hourly/nogaps0.5c-sea_${y}_03hr_TaqaQrQp.D fort.${N} &
    else
      ${pget} ${W}/6hourly/nogaps1.0a-sea_${y}_06hr_TaqaQrQp.D fort.${N} &
    endif
#   ${pget} ${W}/era40-sea_${y}_06hr_TaqaQrQp.D fort.${N} &
  endif
end
C
C --- executable
C
/bin/cp /u/home/wallcraf/hycom/ALL/force/src/ap . &
/bin/cp /u/home/wallcraf/hycom/ALL/force/src/aphf_meanfit . &
wait
chmod ug+rx ap
chmod ug+rx aphf_meanfit
ls -laFq
C
C --- NAMELIST input.
C
touch   fort.05i
/bin/rm fort.05i
cat <<E-o-D  > fort.05i
 &AFTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = 'NOGAPS 0.5c, 3-hrly, ocn-only, MKS',
 /
 &AFTIME
  HMKS   =   1.0,          !kg/kg             to kg/kg
  RMKS   =   1.0,          !W/m**2 into ocean to W/m**2 into ocean
  PMKS   =   1.0,          !m/s    into ocean to m/s    into ocean
  BIASPC =   0.0,
  PCMEAN =   0.0,
  BIASRD =   0.0,
  RDMEAN =   0.0,
  FSTART = ${TS},
  TSTART = ${TS},
  TMAX   = ${TM},
 /
 &AFFLAG
  IFFILE =   5,  !3:monthly-climo; 5:actual-day;
  IFTYPE =   4,  !5:Ta-Ha-Qr-Qp-Pc; 4:Ta-Ha-Qr-Qp; 2:Qr; 1:Pc;
  INTERP =   0,  !0:piecewise-linear; 1:cubic-spline;
  INTMSK =   0,  !0:no mask; 1:land/sea=0/1; 2:land/sea=1/0;
 /
E-o-D
C
C --- run the flux interpolation.
C
touch fort.10 fort.10a
/bin/rm -f fort.1[0-4] fort.1[0-4]a
C
setenv FOR010A fort.10a
setenv FOR011A fort.11a
setenv FOR012A fort.12a
setenv FOR013A fort.13a
setenv FOR014A fort.14a
C
switch ($OS)
case 'SunOS':
case 'Linux':
case 'OSF1':
case 'AIX':
    /bin/rm -f core
    touch core
    ./ap < fort.05i
    breaksw
case 'IRIX64':
    /bin/rm -f core
    touch core
    assign -V
    ./ap < fort.05i
    assign -V
    assign -R
    breaksw
case 'unicos':
    /bin/rm -f core
    touch core
    assign -V
    ./ap < fort.05i
    if (! -z core)  debug -s ap core
    assign -V
    assign -R
    breaksw
endsw
C
C --- Output, use .A and .B if further correction is needed.
C
/bin/mv fort.10  airtmp_${Y01}${A}.b
/bin/mv fort.10a airtmp_${Y01}${A}.a
/bin/mv fort.11  vapmix_${Y01}${A}.b
/bin/mv fort.11a vapmix_${Y01}${A}.a
/bin/mv fort.12  radflx_${Y01}${A}.b
/bin/mv fort.12a radflx_${Y01}${A}.a
/bin/mv fort.13  shwflx_${Y01}${A}.b
/bin/mv fort.13a shwflx_${Y01}${A}.a
if (-e fort.14) then
  /bin/mv fort.14  precip_${Y01}${A}.b
  /bin/mv fort.14a precip_${Y01}${A}.a
endif
C
C --- Annual Mean airtmp correction.
C
if (-e airtmp_${Y01}${A}.A) then
  setenv FOR010A fort.10A
  setenv FOR020A fort.20A
  setenv FOR030A fort.30A
C
  setenv Q era40-fnmoc_clip
  /bin/rm -f fort.10 fort.10A
  /bin/rm -f fort.20 fort.20A fort.30A
  /bin/mv airtmp_${Y01}${A}.B fort.20
  /bin/mv airtmp_${Y01}${A}.A fort.20A
  touch  airtmp_ann_${Q}.a
  if (-z airtmp_ann_${Q}.a) then
    ${pget} ${D}/../../force/offset/airtmp_ann_${Q}.a airtmp_ann_${Q}.a
  endif
  ln -sf airtmp_ann_${Q}.a fort.30A
  ./aphf_meanfit <<'E-o-D'
 &AFMEAN
  TITLE = '1234567890123456789012345678901234567890123456789012345678901234567890123456789',
  TITLE = 'Offset to ERA40 annual mean Ta (-2 < offset < +2)',
  S1     = 1.0,  !no bias
 /
'E-o-D'
  mv fort.10  airtmp_${Y01}${A}.b
  mv fort.10A airtmp_${Y01}${A}.a
  /bin/rm -f fort.20 fort.20A fort.30A
endif
C
C --- Annual Mean Radiative Correction
C
if (-e radflx_${Y01}${A}.A) then
  setenv FOR010A fort.10A
  setenv FOR011A fort.11A
  setenv FOR020A fort.20A
  setenv FOR021A fort.21A
  setenv FOR030A fort.30A
  setenv FOR031A fort.31A
C
  /bin/rm -f fort.10 fort.10A
  /bin/rm -f fort.20 fort.20A fort.30A
  /bin/rm -f fort.21 fort.21A fort.31A
  /bin/mv radflx_${Y01}${A}.B fort.20
  /bin/mv radflx_${Y01}${A}.A fort.20A
  /bin/mv shwflx_${Y01}${A}.B fort.21
  /bin/mv shwflx_${Y01}${A}.A fort.21A
C
  setenv Q fdQfnmoc
  touch  lwflux_ann_${Q}.a
  if (-z lwflux_ann_${Q}.a) then
    ${pget} ${D}/../../force/offset/lwflux_ann_${Q}.a lwflux_ann_${Q}.a &
  endif
  touch  shwflx_ann_${Q}.a
  if (-z shwflx_ann_${Q}.a) then
    ${pget} ${D}/../../force/offset/lwflux_ann_${Q}.a shwflx_ann_${Q}.a &
  endif
  wait
  ln -sf lwflux_ann_${Q}.a fort.30A
  ln -sf shwflx_ann_${Q}.a fort.31A
  ./aphf_meanfit <<'E-o-D'
 &AFMEAN
  TITLE = '1234567890123456789012345678901234567890123456789012345678901234567890123456789',
  TITLE = 'Biased to ISSCP FD annual mean Qsw and Qlw',
          'Biased to ISSCP FD annual mean Qsw',
  SOLAR  = .TRUE.,
  S0     = 0.0,  !no offset for Qlw
           0.0,  !no offset for Qsw
 /
'E-o-D'
  mv fort.10  radflx_${Y01}${A}.b
  mv fort.10A radflx_${Y01}${A}.a
  mv fort.11  shwflx_${Y01}${A}.b
  mv fort.11A shwflx_${Y01}${A}.a
  /bin/rm -f fort.20 fort.20A fort.30A
  /bin/rm -f fort.21 fort.21A fort.31A
endif
C
if (-e ./SAVE) then
  ln airtmp_${Y01}${A}.a ./SAVE/airtmp_${Y01}${A}.a
  ln airtmp_${Y01}${A}.b ./SAVE/airtmp_${Y01}${A}.b
  ln vapmix_${Y01}${A}.a ./SAVE/vapmix_${Y01}${A}.a
  ln vapmix_${Y01}${A}.b ./SAVE/vapmix_${Y01}${A}.b
  ln radflx_${Y01}${A}.a ./SAVE/radflx_${Y01}${A}.a
  ln radflx_${Y01}${A}.b ./SAVE/radflx_${Y01}${A}.b
  ln shwflx_${Y01}${A}.a ./SAVE/shwflx_${Y01}${A}.a
  ln shwflx_${Y01}${A}.b ./SAVE/shwflx_${Y01}${A}.b
  if (-e fort.14) then
    ln precip_${Y01}${A}.a ./SAVE/precip_${Y01}${A}.a
    ln precip_${Y01}${A}.b ./SAVE/precip_${Y01}${A}.b
  endif
endif
C
C --- CICE
C
if (-e ../cice) then
  cd ../cice
  setenv IDM 168
  setenv JDM 280
  /bin/rm vapmix_${Y01}${A}.r airtmp_${Y01}${A}.a airtmp_${Y01}${A}.r netrad_${Y01}${A}.r netQlw_${Y01}${A}.a
  hycom2raw8 ../flux/vapmix_${Y01}${A}.a ${IDM} ${JDM} vapmix_${Y01}${A}.r >! vapmix_${Y01}${A}.B
# airtmp is in degK
  hycom_expr ../flux/airtmp_${Y01}${A}.a ONE ${IDM} ${JDM} 1.0 273.15 airtmp_${Y01}${A}.a >! airtmp_${Y01}${A}.b
  hycom2raw8         airtmp_${Y01}${A}.a     ${IDM} ${JDM}            airtmp_${Y01}${A}.r >! airtmp_${Y01}${A}.B
# netrad may not be needed, create it anyway
  hycom2raw8 ../flux/shwflx_${Y01}${A}.a ${IDM} ${JDM} netrad_${Y01}${A}.r >! netrad_${Y01}${A}.B
# incoming Qlw from netQlw and Tsur, just calculate Qlw here
  hycom_expr ../flux/radflx_${Y01}${A}.a ../flux/shwflx_${Y01}${A}.a ${IDM} ${JDM} 1.0 -1.0 netQlw_${Y01}${A}.a >! netQlw_${Y01}${A}.b
  if (-e ./SAVE) then
    foreach f ( vapmix_${Y01}${A}.[rB] airtmp_${Y01}${A}.[rB] )
      ln ${f} ./SAVE/${f}
    end
  endif
endif
C
C  --- END OF JUST IN TIME FLUX GENERATION SCRIPT.
C
#/usr/lpp/LoadL/full/bin/llq -w $LOADL_STEP_ID
