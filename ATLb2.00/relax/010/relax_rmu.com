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
case 'unicos':
#   newacct NO2030
#   source ~/.login
    breaksw
default:
    echo 'Unknown Operating System: ' $OS
    exit (1)
endsw
C
C --- Generate a HYCOM relaxation mask.
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
case 'Linux':
    setenv pget cp
    setenv pput cp
    breaksw
case 'unicos':
    setenv pget ~wallcraf/bin/pget
    setenv pput ~wallcraf/bin/pput
    breaksw
endsw
C
C --- E is experiment number, from EXPT.src
C --- R is region identifier, from EXPT.src
C --- T is topog. identifier, from EXPT.src
C
C --- P is primary path,
C --- C is scratch directory,
C --- D is permanent directory,
C
source EXPT.src
C
setenv P hycom/${R}/relax/${E}
C
switch ($OS)
case 'SunOS':
    if (-e /scr) then
      setenv C /scr/${user}/$P
      setenv D            ~/$P
    else
      setenv C /net/hermes/scrb/${user}/$P
      setenv D                        ~/$P
    endif
    breaksw
case 'Linux':
    if (-e /external/fast) then
      setenv C /external/fast/${user}/$P
      setenv D                      ~/$P
    else
      setenv C ~/$P/SCRATCH
      setenv D ~/$P
    endif
    breaksw
case 'unicos':
case 'unicos-t3d':
    setenv C /tmp/${user}/$P
    setenv D            ~/$P
    breaksw
endsw
C
mkdir -p $C
cd       $C
C
touch   rmu fort.51 fort.51A regional.grid.a regional.grid.b
/bin/rm rmu fort.51 fort.51A regional.grid.a regional.grid.b
C
C --- Input.
C
touch fort.51 fort.51A
if (-z fort.51) then
  ${pget} ${D}/../../topo/depth_${R}_${T}.b fort.51 &
endif
if (-z fort.51A) then
  ${pget} ${D}/../../topo/depth_${R}_${T}.a fort.51A &
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
touch rmu
if (-z rmu) then
# ${pget} ${D}/../../../ALL/relax/src/rmu . &
  cp      ${D}/../../../ALL/relax/src/rmu . &
endif
wait
chmod a+rx rmu
C
switch ($OS)
case 'SunOS':
    /bin/rm -f core rmu.trace
    touch core
case 'Linux':
    /bin/rm -f core
    touch core
case 'unicos':
    /bin/rm -f core
    touch core
    setenv NCPUS 1
endsw
/bin/rm -f fort.21 fort.21A
C
setenv FOR021A fort.21A
setenv FOR051A fort.51A
C
./rmu <<'E-o-D'
 &MASK
  CTITLE = 'N and S boundary: 5 grid points with 20 to 120 day e-folding time',
  IF     =   1,     1,     1,     1,     1,
             1,     1,     1,     1,     1,
  IL     =  56,    56,    56,    56,    56,
            56,    56,    56,    56,    56,
  JF     =   1,     2,     3,     4,     5,
            51,    50,    49,    48,    47,
  JL     =   1,     2,     3,     4,     5,
            51,    50,    49,    48,    47,
  EFOLD  =  20.0,  45.0,  70.0,  95.0, 120.0,
            20.0,  45.0,  70.0,  95.0, 120.0,
 /
'E-o-D'
C
switch ($OS)
case 'SunOS':
    if (-e rmu.trace) then
      cat rmu.trace
    endif
    breaksw
case 'Linux':
    breaksw
case 'unicos':
    if (! -z core) then
      debug -s rmu core
    endif
    breaksw
endsw
C
C --- Output.
C
${pput} fort.21  ${D}/relax_rmu.b
${pput} fort.21A ${D}/relax_rmu.a
#C
#C --- Delete all scratch directory files.
#C
#/bin/rm -f *
