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
C --- create a HYCOM precip_zero file.
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
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
endif
C
#${pget} ${D}/../../../ALL/force/src/pcip_zero . &
cp ${D}/../../../ALL/force/src/pcip_zero . &
wait
chmod a+rx pcip_zero
C
/bin/rm -f fort.13 fort.13A
C
setenv FOR013A fort.13A
C
switch ($OS)
case 'SunOS':
    /bin/rm -f core pcip_zero.trace
    touch core
default:
    /bin/rm -f core
    touch core
endsw
C
./pcip_zero
C
switch ($OS)
case 'SunOS':
    if (-e pcip_zero.trace) cat pcip_zero.trace
    breaksw
default:
endsw
C
C --- Output.
C
${pput} fort.13  ${D}/precip_zero.b
${pput} fort.13A ${D}/precip_zero.a
C
C --- Delete all files.
C
#/bin/rm -f *
