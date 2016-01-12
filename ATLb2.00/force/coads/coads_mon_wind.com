#
#@$-lt 0:04:00
#@$-lT 0:04:00
#@$-lm 4mw
#@$-lM 4mw
#@$-s  /bin/csh
#@$
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
case 'unicosmk':
    breaksw
default:
    echo 'Unknown Operating System: ' $OS
    exit (1)
endsw
C
C --- Create a wind stress data file.
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
case 'XXX':
    setenv pget ~wallcraf/bin/pget
    setenv pput ~wallcraf/bin/pput
    breaksw
default:
    setenv pget cp
    setenv pput cp
endsw
C
C --- P is primary path,
C --- S is scratch directory,
C --- D is permanent directory,
C
setenv P hycom/ATLb2.00/force/coads
C
switch ($OS)
case 'SunOS':
    if (-e /scr) then
#                  NAVO MSRC
      setenv S /scr/${user}/$P
      setenv D            ~/$P
    else if (-e /net/hermes/scrb) then
#                  NRLSSC
      setenv S /net/hermes/scrb/${user}/$P
      setenv D                        ~/$P
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
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
      setenv D ~/$P
    endif
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
touch      fort.44 fort.44A fort.45 fort.45A fort.71
/bin/rm -f fort.44 fort.44A fort.45 fort.45A fort.71
${pget}                       $D/../offset/tauewd_zero.b fort.44  &
${pget}                       $D/../offset/tauewd_zero.a fort.44A &
${pget}                       $D/../offset/taunwd_zero.b fort.45  &
${pget}                       $D/../offset/taunwd_zero.a fort.45A &
${pget} ~wallcraf/wind_ieee/coads/uwm_coads_monmn_unsm.d fort.71  &
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
endif
C
#${pget} ${D}/../../../ALL/force/src/wi_100_co . &
cp ${D}/../../../ALL/force/src/wi_100_co . &
wait
chmod a+rx wi_100_co
C
/bin/rm -f fort.1[012]*
C
setenv FOR010A fort.10A
setenv FOR011A fort.11A
setenv FOR012A fort.12A
C
setenv FOR044A fort.44A
setenv FOR045A fort.45A
C
switch ($OS)
case 'SunOS':
    /bin/rm -f core wi_100_co.trace
    touch core
    breaksw
default:
    /bin/rm -f core
    touch core
endsw
/bin/rm -f fort.10 fort.11 fort.12 fort.46 fort.47 fort.48
./wi_100_co <<'E-o-D'
 &WWTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = 'COADS monthly, MKS',
 /
 &WWTIME
  SPDMIN = 0.0,
 /
 &WWFLAG
  IGRID  = 1,
  INTERP = 1,
  IWFILE = 1,
 /
'E-o-D'
switch ($OS)
case 'SunOS':
    if (-e wi_100_co.trace) cat wi_100_co.trace
    breaksw
default:
endsw
C
C --- Output.
C
${pput} fort.10  ${D}/tauewd.b
${pput} fort.10A ${D}/tauewd.a
${pput} fort.11  ${D}/taunwd.b
${pput} fort.11A ${D}/taunwd.a
${pput} fort.12  ${D}/wndspd.b
${pput} fort.12A ${D}/wndspd.a
C
C --- Delete all files.
C
#/bin/rm -f *
