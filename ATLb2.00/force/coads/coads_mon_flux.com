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
C --- Create HYCOM heat flux files.
C --- Zero W/m**2 radiation correction, +ve means gained by ocean.
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
touch      fort.71
/bin/rm -f fort.71
${pget} ~wallcraf/flux_ieee/coads/coads_mon_taqaqrqppc.d fort.71 &
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
endif
C
#${pget} ${D}/../../../ALL/force/src/ap_100_co . &
cp ${D}/../../../ALL/force/src/ap_100_co . &
wait
chmod a+rx ap_100_co
C
/bin/rm -f fort.1[01234]*
C
setenv FOR010A fort.10A
setenv FOR011A fort.11A
setenv FOR012A fort.12A
setenv FOR013A fort.13A
setenv FOR014A fort.14A
C
switch ($OS)
case 'SunOS':
    /bin/rm -f core ap_100_co.trace
    touch core
default:
    /bin/rm -f core
    touch core
endsw
C
./ap_100_co <<'E-o-D'
 &AFTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = 'COADS monthly, MKS',
 /
 &AFTIME
  BIASPC =   0.0,
  BIASRD =   0.0,
 /
 &AFFLAG
  IFTYPE =   5,
  INTERP =   1,
 /
'E-o-D'
C
switch ($OS)
case 'SunOS':
    if (-e ap_100_co.trace) cat ap_100_co.trace
    breaksw
default:
endsw
C
C --- Output.
C
${pput} fort.10  ${D}/airtmp.b
${pput} fort.10A ${D}/airtmp.a
${pput} fort.11  ${D}/vapmix.b
${pput} fort.11A ${D}/vapmix.a
${pput} fort.12  ${D}/radflx.b
${pput} fort.12A ${D}/radflx.a
${pput} fort.13  ${D}/shwflx.b
${pput} fort.13A ${D}/shwflx.a
${pput} fort.14  ${D}/precip.b
${pput} fort.14A ${D}/precip.a
C
C --- Delete all files.
C
#/bin/rm -f *
