#
#@ job_name         = 010_tsur_MN
#@ output           = $(job_name).log
#@ error            = $(job_name).log
#@ restart          = yes
#@ job_type         = serial
#@ node_usage       = shared
#@ wall_clock_limit = 6:00:00
#@ resources        = ConsumableCpus(1) ConsumableMemory(26gb)
#@ account_no       = NRLSS018
#@ class            = bigmem
#@ queue
#
set echo
set time=1
#
# --- extract 2-d fields along temperature surfaces
# --- from a single HYCOM mean archive file.
#
# --- output can be formatted, unformatted (BINARY), .[ab] (HYCOM).
# --- or use archv2ncdf2d for netCDF output.
#
# --- output is HYCOM .a files.
#
setenv R ATLb2.00
setenv T 01
setenv E  010
setenv X 01.0
setenv Y 0020
setenv A ""
#
cd /net/ajax/data/$user/hycom/${R}/expt_${X}/data/meanstd
#cd /scr/$user/hycom/${R}/expt_${X}/data/meanstd
#
touch regional.depth.a regional.depth.b
if (-z regional.depth.a) then
  /bin/rm regional.depth.a
  /bin/ln -s ~/hycom/${R}/topo/depth_${R}_${T}.a regional.depth.a
endif
if (-z regional.depth.b) then
  /bin/rm regional.depth.b
  /bin/ln -s ~/hycom/${R}/topo/depth_${R}_${T}.b regional.depth.b
endif
#
touch regional.grid.a regional.grid.b
if (-z regional.grid.a) then
  /bin/rm regional.grid.a
  /bin/ln -s ~/hycom/${R}/topo/regional.grid.a .
endif
if (-z regional.grid.b) then
  /bin/rm regional.grid.b
  /bin/ln -s ~/hycom/${R}/topo/regional.grid.b .
endif
#
    setenv FOR051A ${E}_archMN.${Y}${A}_bot.a
    setenv FOR051  ${E}_archMN.${Y}${A}_bot.b
    setenv FOR031A ${E}_archMN.${Y}${A}_3tl.a
    setenv FOR031  ${E}_archMN.${Y}${A}_3tl.b
    setenv FOR032A ${E}_archMN.${Y}${A}_3td.a
    setenv FOR032  ${E}_archMN.${Y}${A}_3td.b
    setenv FOR033A ${E}_archMN.${Y}${A}_3tt.a
    setenv FOR033  ${E}_archMN.${Y}${A}_3tt.b
    setenv FOR034A ${E}_archMN.${Y}${A}_3ts.a
    setenv FOR034  ${E}_archMN.${Y}${A}_3ts.b
    setenv FOR035A ${E}_archMN.${Y}${A}_3tr.a
    setenv FOR035  ${E}_archMN.${Y}${A}_3tr.b
    /bin/rm $FOR051
    /bin/rm $FOR051A
    /bin/rm $FOR031  $FOR032  $FOR033  $FOR034  $FOR035
    /bin/rm $FOR031A $FOR032A $FOR033A $FOR034A $FOR035A
    ~wallcraf/hycom/ALL/archive/src/archv2data2t <<E-o-D
archMN.${Y}${A}.a
HYCOM
 000	'iexpt ' = experiment number x10 (000=from archive file)
   0	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
 57	'idm   ' = longitudinal array size
 52	'jdm   ' = latitudinal  array size
   0	'itest ' = longitudinal test point (optional, default 0)
   0	'jtest ' = latitudinal  test point (optional, default 0)
 22	'kdm   ' = number of layers
 25.0	'thbase' = reference density (sigma units)
   0	'smooth' = smooth the layered fields (0=F,1=T)
   1	'iorign' = i-origin of plotted subregion
   1	'jorign' = j-origin of plotted subregion
   0	'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
   0	'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
  13	'ktemp ' = number of temperature surfaces to sample
  24.5	'tsur  ' = sample temperaure
  19.5	'tsur  ' = sample temperaure
  17.0	'tsur  ' = sample temperaure
  14.5	'tsur  ' = sample temperaure
   7.0	'tsur  ' = sample temperaure
   4.0	'tsur  ' = sample temperaure
   3.6	'tsur  ' = sample temperaure
   2.8	'tsur  ' = sample temperaure
   2.4	'tsur  ' = sample temperaure
   2.0	'tsur  ' = sample temperaure
   1.6	'tsur  ' = sample temperaure
   1.3	'tsur  ' = sample temperaure
   1.0	'tsur  ' = sample temperaure
   0	'botio ' = bathymetry       I/O unit (0 no I/O)
  31	'layio ' = layer k   i.dep. I/O unit (0 no I/O)
  32	'depio ' = layer k   thick. I/O unit (0 no I/O)
  33	'temio ' = layer k   temp   I/O unit (0 no I/O)
  34	'salio ' = layer k   saln.  I/O unit (0 no I/O)
  35	'tthio ' = layer k   dens,  I/O unit (0 no I/O)
E-o-D
#
#/usr/lpp/LoadL/full/bin/llq -w $LOADL_STEP_ID
