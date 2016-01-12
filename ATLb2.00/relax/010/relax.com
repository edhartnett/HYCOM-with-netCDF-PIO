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
C --- Convert z-level climatology to HYCOM layers.
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
touch   relax_sig0 blkdat.input fort.51 fort.51A
/bin/rm relax_sig0 blkdat.input fort.51 fort.51A
C
C --- 12 months
C
foreach MM ( 01 02 03 04 05 06 07 08 09 10 11 12 )
C
C --- Input.
C
touch      fort.71 fort.71A fort.72 fort.72A
/bin/rm -f fort.71 fort.71A fort.72 fort.72A
${pget} ${D}/../levitus/dens_sig0_m${MM}.b fort.71  &
${pget} ${D}/../levitus/dens_sig0_m${MM}.a fort.71A &
${pget} ${D}/../levitus/temp_sig0_m${MM}.b fort.72  &
${pget} ${D}/../levitus/temp_sig0_m${MM}.a fort.72A &
touch fort.51 fort.51A
if (-z fort.51) then
  ${pget} ${D}/../../topo/depth_${R}_${T}.b fort.51  &
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
touch blkdat.input
if (-z blkdat.input) then
  ${pget} ${D}/blkdat.input blkdat.input &
endif
C
touch relax_sig0
if (-z relax_sig0) then
# ${pget} ${D}/../../../ALL/relax/src/relax_sig0 . &
  cp      ${D}/../../../ALL/relax/src/relax_sig0 . &
endif
wait
chmod a+rx relax_sig0
C
sed -e "s/^[ 	0-9]*'month ' =/  ${MM}	  'month ' =/" blkdat.input >! fort.99
C
switch ($OS)
case 'SunOS':
    /bin/rm -f core relax_sig0.trace
    touch core
    breaksw
default:
    /bin/rm -f core
    touch core
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
./relax_sig0
switch ($OS)
case 'SunOS':
    if (-e relax_sig0.trace) then
      cat relax_sig0.trace
    endif
    breaksw
default:
endsw
C
C --- Output.
C
mv fort.10  relax_tem_m${MM}.b
mv fort.10A relax_tem_m${MM}.a
mv fort.11  relax_sal_m${MM}.b
mv fort.11A relax_sal_m${MM}.a
mv fort.12  relax_int_m${MM}.b
mv fort.12A relax_int_m${MM}.a
C
setenv DAYM `echo ${MM} | awk '{printf("0000_%3.3d_00\n",30*($1-1)+16)}'`
#setenv DAYM `echo ${MM} | awk '{printf("0000_%3.3d_00\n",30.5*($1-1)+16)}'`
#setenv DAYM `echo ${MM} | awk '{printf("0000_%3.3d_00\n",30.5*($1-1)+1)}'`
${pput} fort.21  ${D}/relax.${DAYM}.b
${pput} fort.21A ${D}/relax.${DAYM}.a
C
C --- end of month foreach loop
C
/bin/rm fort.7[12]
end
#
# --- Merge monthly climatologies into one file.
#
cp relax_int_m01.b relax_int.b
cp relax_sal_m01.b relax_sal.b
cp relax_tem_m01.b relax_tem.b
#
foreach MM ( 02 03 04 05 06 07 08 09 10 11 12 )
  tail +6 relax_int_m${MM}.b >> relax_int.b
  tail +6 relax_sal_m${MM}.b >> relax_sal.b
  tail +6 relax_tem_m${MM}.b >> relax_tem.b
end
#
cp relax_int_m01.a relax_int.a
cp relax_sal_m01.a relax_sal.a
cp relax_tem_m01.a relax_tem.a
#
foreach MM ( 02 03 04 05 06 07 08 09 10 11 12 )
  cat relax_int_m${MM}.a >> relax_int.a
  cat relax_sal_m${MM}.a >> relax_sal.a
  cat relax_tem_m${MM}.a >> relax_tem.a
end
${pput} relax_int.b ${D}/relax_int.b
${pput} relax_int.a ${D}/relax_int.a
${pput} relax_sal.b ${D}/relax_sal.b
${pput} relax_sal.a ${D}/relax_sal.a
${pput} relax_tem.b ${D}/relax_tem.b
${pput} relax_tem.a ${D}/relax_tem.a
#
# --- delete the monthly files
#
/bin/rm relax_int_m??.[ab]
/bin/rm relax_sal_m??.[ab]
/bin/rm relax_tem_m??.[ab]
#C
#C --- Delete all scratch directory files.
#C
#/bin/rm -f *
