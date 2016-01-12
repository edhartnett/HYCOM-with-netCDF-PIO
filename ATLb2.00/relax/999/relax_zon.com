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
case 'unicos':
#   newacct NO2030
#   source ~/.login
    breaksw
default:
    echo 'Unknown Operating System: ' $OS
    exit (1)
endsw
C
C --- Summer zonal climatological statistics.
C
switch ($OS)
case 'SunOS':
    setenv pget cp
    setenv pput cp
    breaksw
case 'unicos':
    setenv pget ~wallcraf/bin/pget
    setenv pput ~wallcraf/bin/pput
    breaksw
endsw
C
C --- E is experiment number
C --- P is primary path,
C --- C is scratch directory,
C --- D is permanent directory,
C
setenv E 999
setenv P hycom/ATLa2.00/relax/${E}
C
switch ($OS)
case 'SunOS':
    setenv C /net/hermes/scrb/${user}/$P
    setenv D                        ~/$P
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
touch   relax_zon blkdat.input fort.51
/bin/rm relax_zon blkdat.input fort.51
C
C --- July.
C
foreach MM ( 07 )
C
C --- Input.
C
touch      fort.71 fort.72
/bin/rm -f fort.71 fort.72
${pget} ${D}/../levitus/dens_m${MM} fort.71 &
${pget} ${D}/../levitus/temp_m${MM} fort.72 &
C
touch fort.51
if (-z fort.51) then
  ${pget} ${D}/../../topo/depth_ATLa2.00_01 fort.51 &
endif
C
touch blkdat.input
if (-z blkdat.input) then
  ${pget} ${D}/blkdat.input blkdat.input &
endif
C
touch relax_zon
if (-z relax_zon) then
  ${pget} ${D}/../src/relax_zon . &
endif
wait
chmod a+rx relax_zon
C
sed -e "s/^[ 	0-9]*'month ' =/  ${MM}	  'month ' =/" blkdat.input >! fort.99
C
switch ($OS)
case 'SunOS':
    /bin/rm -f core relax_zon.trace
    touch core
case 'unicos':
    /bin/rm -f core
    touch core
    setenv NCPUS 1
endsw
/bin/rm -f fort.10 fort.11 fort.12 fort.21
./relax_zon
switch ($OS)
case 'SunOS':
    if (-e relax_zon.trace) then
      cat relax_zon.trace
    endif
    breaksw
case 'unicos':
    ./relax_zon
    if (! -z core) then
      debug -s relax_zon core
    endif
    breaksw
endsw
C
C --- Output.
C
#${pput} fort.10 ${D}/relax_tem_m${MM}
#${pput} fort.11 ${D}/relax_sal_m${MM}
#${pput} fort.12 ${D}/relax_int_m${MM}
C
setenv DAYM `echo ${MM} | awk '{printf("%6.6d\n",30*($1-1))}'`
#setenv DAYM `echo ${MM} | awk '{printf("%6.6d\n",30.5*($1-1))}'`
#setenv DAYM `echo ${MM} | awk '{printf("%6.6d\n",30.5*($1-1)+15)}'`
${pput} fort.21 ${D}/relax.${DAYM}
C
C --- end of month foreach loop
C
/bin/rm fort.7[12]
end
#C
#C --- Delete all scratch directory files.
#C
#/bin/rm -f *
