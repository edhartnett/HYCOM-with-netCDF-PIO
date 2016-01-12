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
C --- Create a monthly climatology wind stress data file.
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
touch      fort.44 fort.44A fort.45 fort.45A fort.71
/bin/rm -f fort.44 fort.44A fort.45 fort.45A fort.71
${pget}              $D/../offset/tauewd_zero.b fort.44  &
${pget}              $D/../offset/tauewd_zero.a fort.44A &
${pget}              $D/../offset/taunwd_zero.b fort.45  &
${pget}              $D/../offset/taunwd_zero.a fort.45A &
ln -s   /net/hermes/data5/metzger/wind/${N}/${W}_strblk.D fort.71  &
#ln -s   ~metzger/force/${N}/monthly/${W}_strblk.D fort.71  &
#${pget} ~metzger/force/${N}/monthly/${W}_strblk.D fort.71  &
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
endif
C
#${pget} ${D}/../../../ALL/force/src/wi . &
cp ${D}/../../../ALL/force/src/wi . &
wait
chmod a+rx wi
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
    /bin/rm -f core wi.trace
    touch core
    breaksw
default:
    /bin/rm -f core
    touch core
endsw
/bin/rm -f fort.10 fort.11 fort.12 fort.46 fort.47 fort.48
./wi <<E-o-D
 &WWTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = "${W} strblk",
 /
 &WWTIME
  WSCALE = 1.0,  !scale factor to mks
  SPDMIN = 0.0,  !minimum allowed wind speed
 /
 &WWFLAG
  IGRID  = 2,  !1:u&v; 2:p-grid
  INTERP = 1,  !0:bilinear; 1:cubic spline
  IWFILE = 1,  !1:clim
  ISPEED = 0,  !0:none; 1:const; 2:kara; 3:coare
  INTMSK = 0,  !0:no mask; 1:land/sea=0/1; 2:land/sea=1/0;
 /
E-o-D
switch ($OS)
case 'SunOS':
    if (-e wi.trace) cat wi.trace
    breaksw
default:
endsw
C
C --- Output.
C
mv fort.10  tauewd.b
mv fort.10A tauewd.a
mv fort.11  taunwd.b
mv fort.11A taunwd.a
C
foreach f ( tauewd.? taunwd.? )
  ${pput} ${f} ${D}/${f} &
end
wait
C
C --- Delete all files.
C
/bin/rm -f fort.* wi
