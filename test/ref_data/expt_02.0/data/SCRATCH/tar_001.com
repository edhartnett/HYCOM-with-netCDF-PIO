#!/bin/csh
#
set echo
set time = 1
set timestamp
C
C --- tar archive files in a batch job
C
setenv OS `uname`
C
C --- pget, pput "copy" files between scratch and permanent storage.
C --- Can both be cp if the permanent filesystem is mounted locally.
C
switch ($OS)
case 'SunOS':
case 'Linux':
#case 'IRIX64':
case 'OSF1':
case 'AIX':
case 'unicos':
    if      (-e ~wallcraf/bin/pget_NAVO) then
      setenv pget ~wallcraf/bin/pget_NAVO
      setenv pput ~wallcraf/bin/pput_NAVO
    else if (-e ~wallcraf/bin/pget) then
      setenv pget ~wallcraf/bin/pget
      setenv pput ~wallcraf/bin/pput
    else
      setenv pget /bin/cp
      setenv pput /bin/cp
    endif
    breaksw
default:
    setenv pget /bin/cp
    setenv pput /bin/cp
endsw
C
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C --- GTAR location of gtar, or any tar supporting "-z" for gzip
C
setenv E 020
setenv P hycom/ATLb2.00/expt_02.0/data
setenv D ~/$P
C
switch ($OS)
case 'HP-UX':
    setenv S     /scratch3/${user}/$P
    breaksw
case 'OSF1':
#                ASC MSRC
    setenv S     /workspace/${user}/$P
    breaksw
case 'SunOS':
    if (-e /net/hermes/scrb) then
#                  NRLSSC
      setenv S     /net/hermes/scrb/${user}/$P
    else
#                  NAVO MSRC
      setenv S     /scr/${user}/$P
    endif
    breaksw
case 'IRIX64':
    setenv S     /scr/${user}/$P
    breaksw
case 'AIX':
    if      (-e /gpfs/work) then
#                  ERDC MSRC, under PBS
      setenv S     /gpfs/work/${user}/$P
      setenv POE  pbspoe
    else if (-e /scr) then
#                  NAVO MSRC, under LoadLeveler
      setenv S     /scr/${user}/$P
      setenv POE  poe
    else if (-e /scratch/tempest) then
#                  MHPCC DC, under LoadLeveler
      setenv S     /scratch/tempest/users/${user}/$P
      setenv POE  poe
    else
#                  ARL MSRC, under GRD
      setenv S     /usr/var/tmp/${user}/$P
      setenv POE  grd_poe
    endif
    breaksw
case 'unicos':
case 'Linux':
    setenv GTAR /bin/tar
    if (-e /work) then
#                  ERDC
      setenv S /work/${user}/$P
    else if (-e /scr) then
#                  NAVO
      setenv S /scr/${user}/$P
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
    endif
    breaksw
endsw
C
setenv A ""
setenv Y01 "001"
C
C --- run in the tar directories.
C
foreach t ( v m e )
  chmod 750 $S/tar${t}_${Y01}${A}
  cd        $S/tar${t}_${Y01}${A}
C
  /bin/rm  -f  ../${E}_arch${t}_${Y01}${A}.tar ${E}_arch?.dummy.*
  if ( $?GTAR ) then
    $GTAR czvf ../${E}_arch${t}_${Y01}${A}.tar.gz .
  else
    /bin/tar cvf ../${E}_arch${t}_${Y01}${A}.tar .
    gzip         ../${E}_arch${t}_${Y01}${A}.tar
  endif
C
C --- clean up in the primary scratch directory.
C
  cd ..
  /bin/rm -r tar${t}_${Y01}${A}
  ${pput} ${E}_arch${t}_${Y01}${A}.tar.gz $D/${E}_arch${t}_${Y01}${A}.tar.gz
end
