#! /bin/csh
#BSUB -P NRLSS018
#BSUB -n 1
#BSUB -W 8:00
#
set echo
set time = 1
C
setenv OS `uname`
switch ($OS)
case 'SunOS':
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
#   newacct NO2030
#   source ~/.login
    breaksw
default:
    echo 'Unknown Operating System: ' $OS
    exit (1)
endsw
C
C --- Create HYCOM precip data file from a monthly atmos. climatology.
C
switch ($OS)
case 'SunOS':
    if (-e /scr) then
      setenv pget ~wallcraf/bin/pget
      setenv pput ~wallcraf/bin/pput
    else
      setenv pget cp
      setenv pput cp
    endif
    breaksw
    breaksw
case 'AIX':
    setenv pget ~wallcraf/bin/pget
    setenv pput ~wallcraf/bin/pput
    breaksw
default:
    setenv pget cp
    setenv pput cp
endsw
C
C --- R is region,            e.g. GLBgx1v3
C --- W is wind dataset name, e.g. era40-sea_1978-2002_mon
C --- N is wind dataset,      e.g. era40
C
setenv R `echo $cwd | sed -e "s?^.*hycom/??" | awk -F"/" '{print $1}'`
setenv W `echo $cwd                          | awk -F"/" '{print $NF}'`
setenv N `echo $W                            | awk -F"-" '{print $1}'`
C
C --- P is primary path,
C --- S is scratch directory,
C --- D is permanent directory,
C
setenv P hycom/${R}/force/$W
C
switch ($OS)
case 'SunOS':
    if (-e /scr) then
#                  NAVO MSRC
      setenv S /scr/${user}/$P
      setenv D            ~/$P
    else if (-e /net/ajax/data) then
#                  NRLSSC
      setenv S /net/ajax/data/${user}/$P
      setenv D                      ~/$P
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
      setenv D ~/$P
    endif
    breaksw
case 'Linux':
    if (-e /external/fast) then
#                  NRLSSC
      setenv S /external/fast/${user}/$P
      setenv D                      ~/$P
    else if (-e /net/ajax/data) then
#                  NRLSSC
      setenv S /net/ajax/data/${user}/$P
      setenv D                      ~/$P
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
      setenv D ~/$P
    endif
    breaksw
case 'AIX':
    if      (-e /gpfs/work) then
#                  ERDC MSRC
      mkdir        /gpfs/work/${user}
      chmod a+rx   /gpfs/work/${user}
      setenv S     /gpfs/work/${user}/$P
      setenv D                      ~/$P
    else if (-e /scr) then
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv D                ~/$P
    else
#                  ARL MSRC
      mkdir        /usr/var/tmp/${user}
      chmod a+rx   /usr/var/tmp/${user}
      setenv S     /usr/var/tmp/${user}/$P
      setenv D                        ~/$P
    endif
    breaksw
case 'unicos':
case 'unicos-t3d':
    setenv C /tmp/${user}/$P
    setenv S /tmp/${user}/$P
    setenv D            ~/$P
    setenv F            ~/$P
    breaksw
default:
#                  Single Disk
    setenv S ~/$P/SCRATCH
    setenv D ~/$P
endsw
C
mkdir -p $S
cd       $S
C
C --- Input.
C
touch      fort.71
/bin/rm -f fort.71
ln -s    ~/force/${N}/monthly/${W}_ttlpcp.D fort.71  &
#${pget} ~metzger/force/${N}/monthly/${W}_ttlpcp.D fort.71  &
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  cp ${D}/../../topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  cp ${D}/../../topo/regional.grid.a regional.grid.a &
endif
C
#${pget} ${D}/../../../ALL/force/src/kp . &
cp ${D}/../../../ALL/force/src/kp . &
C
wait
chmod a+rx kp
C
/bin/rm -f fort.1[01234]*
C
setenv FOR010A fort.10A
C
switch ($OS)
case 'SunOS':
    /bin/rm -f core kp.trace
    touch core
case 'unicos':
    /bin/rm -f core
    touch core
    setenv NCPUS 1
endsw
/bin/rm -f fort.10
./kp <<E-o-D
 &AFTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = "${W} ttlpcp",
  CNAME  = 'precip',
 &END
 &AFTIME
  PARMIN =    0.0,  !precip is positive into ocean
  PARMAX =   99.0,  !disable maximum
  PAROFF =    0.0,  !disable offset
  PARSCL =    1.0,  !disable scale factor, already m/s into ocean (ERA40)
  PARSCL =    1.1574074E-5,   !m/day into ocean to m/s into ocean (ERA15)
 &END
 &AFFLAG
  INTERP =   0,  !0:piecewise-linear; 1:cubic-spline;
  INTMSK =   0,  !0:no mask; 1:land/sea=0/1; 2:land/sea=1/0;
 &END
E-o-D
switch ($OS)
case 'SunOS':
    if (-e kp.trace) cat kp.trace
    breaksw
case 'unicos':
    ./kp
    if (! -z core) debug -s kp core
    breaksw
endsw
C
C --- Output.
C
/bin/mv fort.10  ./precip.b
/bin/mv fort.10A ./precip.a
C
${pput} precip.b ${D}/precip.b
${pput} precip.a ${D}/precip.a
C
C --- Delete all files.
C
/bin/rm -f fort.* ./kp
