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
C --- Extract bottom density from z-level climatology.
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
C --- S is scratch directory,
C --- D is permanent directory,
C
source EXPT.src
C
setenv P hycom/${R}/relax/${E}
C
switch ($OS)
case 'SunOS':
    if (-e /scr) then
      setenv S /scr/${user}/$P
      setenv D            ~/$P
    else
      setenv S /net/hermes/scrb/${user}/$P
      setenv D                        ~/$P
    endif
    breaksw
case 'unicos':
case 'unicos-t3d':
    setenv S /tmp/${user}/$P
    setenv D            ~/$P
    breaksw
endsw
C
mkdir -p $S
cd       $S
C
touch   bottom_sig0 blkdat.input fort.51 fort.51A
/bin/rm bottom_sig0 blkdat.input fort.51 fort.51A
C
C --- representative month
C
foreach MM ( 01 )
C
C --- Input.
C
touch      fort.71 fort.71A fort.72 fort.72A
/bin/rm -f fort.71 fort.71A fort.72 fort.72A
${pget} ${D}/../levitus/dens_sig0_m${MM}.b fort.71  &
${pget} ${D}/../levitus/dens_sig0_m${MM}.a fort.71A &
${pget} ${D}/../levitus/temp_sig0_m${MM}.b fort.72  &
${pget} ${D}/../levitus/temp_sig0_m${MM}.a fort.72A &
C
touch fort.51 fort.51A
if (-z fort.51) then
  ${pget} ${D}/../../topo/depth_${R}_${T}.b fort.51  &
endif
if (-z fort.51A) then
  ${pget} ${D}/../../topo/depth_${R}_${T}.a fort.51A &
endif
C
touch blkdat.input
if (-z blkdat.input) then
  ${pget} ${D}/blkdat.input blkdat.input &
endif
C
touch bottom_sig0
if (-z bottom_sig0) then
# ${pget} ${D}/../src/bottom_sig0 . &
  cp      ${D}/../src/bottom_sig0 . &
endif
wait
chmod a+rx bottom_sig0
C
sed -e "s/^[ 	0-9]*'month ' =/  ${MM}	  'month ' =/" blkdat.input >! fort.99
C
switch ($OS)
case 'SunOS':
    /bin/rm -f core bottom_sig0.trace
    touch core
case 'unicos':
    /bin/rm -f core
    touch core
    setenv NCPUS 1
endsw
C
setenv FOR010A fort.10A
setenv FOR011A fort.11A
setenv FOR012A fort.12A
setenv FOR021A fort.21A
setenv FOR051A fort.51A
setenv FOR071A fort.71A
setenv FOR072A fort.72A
/bin/rm -f fort.10  fort.11  fort.12  fort.21
/bin/rm -f fort.10A fort.11A fort.12A fort.21A
C
./bottom_sig0
switch ($OS)
case 'SunOS':
    if (-e bottom_sig0.trace) then
      cat bottom_sig0.trace
    endif
    breaksw
case 'unicos':
    ./bottom_sig0
    if (! -z core) then
      debug -s bottom_sig0 core
    endif
    breaksw
endsw
C
C --- Output.
C
mv fort.10  bottom_tem_m${MM}.b
mv fort.10A bottom_tem_m${MM}.a
mv fort.11  bottom_sal_m${MM}.b
mv fort.11A bottom_sal_m${MM}.a
mv fort.12  bottom_int_m${MM}.b
mv fort.12A bottom_int_m${MM}.a
C
setenv DAYM `echo ${MM} | awk '{printf("0000_%3.3d_00\n",30*($1-1)+16)}'`
#setenv DAYM `echo ${MM} | awk '{printf("0000_%3.3d_00\n",30.5*($1-1)+16)}'`
#setenv DAYM `echo ${MM} | awk '{printf("0000_%3.3d_00\n",30.5*($1-1)+1)}'`
${pput} fort.21  ${D}/bottom.${DAYM}.b
${pput} fort.21A ${D}/bottom.${DAYM}.a
C
C --- end of month foreach loop
C
/bin/rm fort.7[12]
end
#
# --- Merge monthly climatologies into one file.
#
cp bottom_int_m01.b bottom_int.b
cp bottom_sal_m01.b bottom_sal.b
cp bottom_tem_m01.b bottom_tem.b
#
#foreach MM ( 02 03 04 05 06 07 08 09 10 11 12 )
#  tail +6 bottom_int_m${MM}.b >> bottom_int.b
#  tail +6 bottom_sal_m${MM}.b >> bottom_sal.b
#  tail +6 bottom_tem_m${MM}.b >> bottom_tem.b
#end
#
cp bottom_int_m01.a bottom_int.a
cp bottom_sal_m01.a bottom_sal.a
cp bottom_tem_m01.a bottom_tem.a
#
#foreach MM ( 02 03 04 05 06 07 08 09 10 11 12 )
#  cat bottom_int_m${MM}.a >> bottom_int.a
#  cat bottom_sal_m${MM}.a >> bottom_sal.a
#  cat bottom_tem_m${MM}.a >> bottom_tem.a
#end
${pput} bottom_int.b ${D}/bottom_int.b
${pput} bottom_int.a ${D}/bottom_int.a
${pput} bottom_sal.b ${D}/bottom_sal.b
${pput} bottom_sal.a ${D}/bottom_sal.a
${pput} bottom_tem.b ${D}/bottom_tem.b
${pput} bottom_tem.a ${D}/bottom_tem.a
#
# --- delete the monthly files
#
/bin/rm bottom_int_m??.[ab]
/bin/rm bottom_sal_m??.[ab]
/bin/rm bottom_tem_m??.[ab]
#C
#C --- Delete all scratch directory files.
#C
#/bin/rm -f *
