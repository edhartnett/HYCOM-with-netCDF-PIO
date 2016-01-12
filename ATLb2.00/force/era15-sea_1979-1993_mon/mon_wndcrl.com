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
C --- Create a wind stress curl data file.
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
      setenv S  /net/ajax/data/${user}/$P
      setenv D                         ~/$P
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
      setenv S  /net/ajax/data/${user}/$P
      setenv D                         ~/$P
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
touch      fort.10 fort.10A fort.11 fort.11A
/bin/rm -f fort.10 fort.10A fort.11 fort.11A
#${pget} $D/tauewd.b fort.10  &
#${pget} $D/tauewd.a fort.10A &
#${pget} $D/taunwd.b fort.11  &
#${pget} $D/taunwd.a fort.11A &
ln -s      tauewd.b fort.10  &
ln -s      tauewd.a fort.10A &
ln -s      taunwd.b fort.11  &
ln -s      taunwd.a fort.11A &
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
endif
C
cp ${D}/../../../ALL/force/src/wi_curl . &
wait
chmod a+rx wi_curl
C
/bin/rm -f fort.12*
C
setenv FOR010A fort.10A
setenv FOR011A fort.11A
setenv FOR012A fort.12A
C
switch ($OS)
case 'SunOS':
    /bin/rm -f core wi_curl.trace
    touch core
    breaksw
default:
    /bin/rm -f core
    touch core
endsw
./wi_curl 
switch ($OS)
case 'SunOS':
    if (-e wi_curl.trace) cat wi_curl.trace
    breaksw
default:
endsw
C
C --- Output.
C
mv fort.12  wndcrl.b
mv fort.12A wndcrl.a
#${pput} wndcrl.b $D/wndcrl.b
#${pput} wndcrl.a $D/wndcrl.a
C
C --- Delete all files.
C
/bin/rm -f wi_curl fort.*
