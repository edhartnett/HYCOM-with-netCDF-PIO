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
C --- Create HYCOM monthly climatology heat flux files.
C --- 0 W/m**2 radiation correction, +ve means gained by ocean.
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
ln -s   ~/force/${N}/monthly/${W}_TaqaQrQp.D fort.71  &
#${pget} ~/force/${N}/monthly/${W}_TaqaQrQp.D fort.71  &
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
endif
C
#${pget} ${D}/../../../ALL/force/src/ap . &
cp ${D}/../../../ALL/force/src/ap . &
wait
chmod a+rx ap
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
    /bin/rm -f core ap.trace
    touch core
default:
    /bin/rm -f core
    touch core
endsw
C
./ap <<E-o-D
 &AFTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = "${W} TaqaQrQp",
 /
 &AFTIME
  BIASPC =   0.0,
  BIASRD =   0.0,
  HMKS   =   1.0,          !kg/kg             to kg/kg
  RMKS   =   1.0,          !W/m**2 into ocean to W/m**2 into ocean
  PMKS   =   1.1574074E-5, !m/day  into ocean to m/s    into ocean
 /
 &AFFLAG
  IFFILE =   3,  !3:monthly-climo; 5:actual-day;
  IFTYPE =   4,  !5:Ta-Ha-Qr-Qp-Pc; 4:Ta-Ha-Qr-Qp; 2:Qr; 1:Pc;
  INTERP =   0,  !0:piecewise-linear; 1:cubic-spline;
  INTMSK =   0,  !0:no mask; 1:land/sea=0/1; 2:land/sea=1/0;
 /
E-o-D
C
switch ($OS)
case 'SunOS':
    if (-e ap.trace) cat ap.trace
    breaksw
default:
endsw
C
C --- Output.
C
/bin/mv fort.10  airtmp.b
/bin/mv fort.10A airtmp.a
/bin/mv fort.11  vapmix.b
/bin/mv fort.11A vapmix.a
/bin/mv fort.12  radflx.b
/bin/mv fort.12A radflx.a
/bin/mv fort.13  shwflx.b
/bin/mv fort.13A shwflx.a
/bin/mv fort.14  precip_zero.b
/bin/mv fort.14A precip_zero.a
C
C
foreach f ( airtmp.? vapmix.? radflx.? shwflx.? precip_zero.? )
  ${pput} ${f} ${D}/${f} &
end
wait
C
C --- Delete all input files.
C
/bin/rm -f fort.* ./ap
