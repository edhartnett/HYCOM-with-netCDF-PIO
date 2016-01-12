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
C --- Interpolate Levitus climtology to a HYCOM region.
C --- KSIGMA=2 for _sig2_ fields
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
C --- D is scalar permanent directory,
C --- L is Levitus directory
C
setenv P hycom/ATLa2.00/relax/levitus
C
switch ($OS)
case 'SunOS':
    if (-e /scr) then
#                  NAVO MSRC
      setenv S /scr/${user}/$P
      setenv D            ~/$P
      setenv L ~/clim/levitus
    else if (-e /net/hermes/scrb) then
#                  NRLSSC
      setenv S /net/hermes/scrb/${user}/$P
      setenv D                        ~/$P
      setenv L ~/clim/levitus
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
      setenv D ~/$P
      setenv L ~/clim/levitus
    endif
    breaksw
case 'Linux':
    if (-e /external/fast) then
#                  NRLSSC
      setenv S /external/fast/${user}/$P
      setenv D                      ~/$P
      setenv L /external/fast/${user}/clim/levitus
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
      setenv D ~/$P
      setenv L ~/clim/levitus
    endif
    breaksw
default:
#                  Single Disk
    setenv S ~/$P/SCRATCH
    setenv D ~/$P
    setenv L ~/clim/levitus
endsw
C
mkdir -p $S
cd       $S
C
touch   z_levitus
/bin/rm z_levitus
C
C --- 12 months
C
foreach MM ( 01 02 03 04 05 06 07 08 09 10 11 12 )
C
C --- Input.
C
touch      fort.71 fort.73
/bin/rm -f fort.71 fort.73
${pget} $L/r_m${MM}.d fort.71 &
${pget} $L/s_m${MM}.d fort.73 &
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
endif
C
touch z_levitus
if (-z z_levitus) then
  ${pget} ${D}/../src/z_levitus . &
endif
wait
chmod a+rx z_levitus
C
switch ($OS)
case 'SunOS':
    /bin/rm -f core z_levitus.trace
    touch core
    breaksw
default:
    /bin/rm -f core
    touch core
endsw
setenv FOR010A fort.10A
setenv FOR011A fort.11A
setenv FOR012A fort.12A
/bin/rm -f fort.10 fort.10A fort.11 fort.11A fort.12 fort.12A
./z_levitus <<E-o-D
 &AFTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = 'Levitus monthly',
 /
 &AFFLAG
  ICTYPE =   2,
  KSIGMA =   2,
  INTERP =   1,
  ITEST  =  40,
  JTEST  =  17,
  MONTH  = $MM,
 /
E-o-D
switch ($OS)
case 'SunOS':
    if (-e z_levitus.trace) then
      cat z_levitus.trace
    endif
    breaksw
default:
endsw
C
C --- Required Output, potential density and temperature.
C
${pput} fort.10  ${D}/temp_sig2_m${MM}.b
${pput} fort.10A ${D}/temp_sig2_m${MM}.a
${pput} fort.12  ${D}/dens_sig2_m${MM}.b
${pput} fort.12A ${D}/dens_sig2_m${MM}.a
C
C --- Optional Output.
C
${pput} fort.11  ${D}/saln_sig2_m${MM}.b
${pput} fort.11A ${D}/saln_sig2_m${MM}.a
C
C --- end of month foreach loop
C
/bin/rm fort.7[123]
end
C
C --- Delete all files.
C
#/bin/rm -f *
