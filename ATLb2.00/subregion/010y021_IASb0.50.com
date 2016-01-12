#!/bin/csh
#
set echo
set time=1
#
# --- form interpolated subregion archive files, ATLb2.00 to IASb0.50.
#
# --- script includes logic to reduce the number of layers,
# --- but it is inactivated here by setting L to "".
#
# --- ALL/bin/hycom_ij2lonlat can be used to find co-located points
# ---  on the two grids.  Since subregion p(1,1) must be on the
# ---  original grid, this is usually the point to reference.
# --- In this case:
# ---   hycom_ij2lonlat 1 1  ~/hycom/IASb0.50/topo/regional.grid.a
# ---    97.000W   3.997N
# ---   hycom_ij2lonlat 1 13 ~/hycom/ATLb2.00/topo/regional.grid.a
# ---    97.000W   3.997N
#
setenv D  ../expt_01.0/data
setenv DT ../topo
setenv R  ../expt_01.0/data/IASb0.50
setenv RT ../../IASb0.50/topo
setenv E  010
setenv Y  021
setenv A  ""
#
# --- single model segment, potentially across two calendar years.
# --- configured for  1-day surface and  6-day 3-d archives.
#
setenv YRFLAG  0
setenv BNSTFQ  1.0
setenv NESTFQ  6.0
#
setenv DF `echo LIMITS | awk -f ${D}/../${E}.awk y01=${Y} ab=${A} | awk '{print $1}'`
setenv DL `echo LIMITS | awk -f ${D}/../${E}.awk y01=${Y} ab=${A} | awk '{print $2}'`
#
touch   tar_list.tmp
/bin/rm tar_list.tmp
touch   tar_list.tmp
#
echo $YRFLAG $BNSTFQ $NESTFQ $DF $DL | ~/hycom/ALL/bin/hycom_nest_dates >! nest_list.tmp
#
foreach ydh ( `paste -d" " -s nest_list.tmp` )
    touch   regional.grid.a
    /bin/rm regional.grid.[ab]
    /bin/ln -s ${DT}/regional.grid.a regional.grid.a
    /bin/ln -s ${DT}/regional.grid.b regional.grid.b
#
    setenv l `cat "${D}/archv.${ydh}.b" | wc -l`
    if (${l} == 32) then
      setenv L ""
    else
#     setenv L "_L22"
      setenv L ""
    endif
    touch   ${R}/archv.${ydh}${L}.b
    /bin/rm ${R}/archv.${ydh}${L}.[ab]
    ~/hycom/ALL/subregion/src/isubregion <<E-o-D
${D}/archv.${ydh}.b
${DT}/depth_ATLb2.00_01.b
${R}/archv.${ydh}${L}.b
${RT}/depth_IASb0.50_03.b
ATLb2.00 interpolated to IASb0.50
  93	  'idm   ' = longitudinal array size
  65	  'jdm   ' = latitudinal  array size
   1	  'irefi ' = longitudinal input  reference location
  13	  'jrefi ' = latitudinal  input  reference location
   1	  'irefo ' = longitudinal output reference location
   1	  'jrefo ' = latitudinal  output reference location
   4	  'ijgrd ' = integer scale factor between input and output grids
   0	  'iceflg' = ice in output archive flag (0=none,1=energy loan model)
   0	  'smooth' = smooth interface depths    (0=F,1=T)
E-o-D
#
    touch    ${R}/archv.${ydh}${L}.b
    if (-z ${R}/archv.${ydh}${L}.b) then
      echo "missing archive file: " ${R}/archv.${ydh}${L}.b
      exit (2)
    endif
#
    if (${L} == "_L22") then
#
# --- reduce any resulting 22-layer archive files (*_L22.[ab]) to
# --- 21-layers using ALL/archive/src/trim_archv.
#
      touch   regional.grid.a    regional.depth.a
      /bin/rm regional.grid.[ab] regional.depth.[ab]
      /bin/ln -s ${RT}/regional.grid.a     regional.grid.a
      /bin/ln -s ${RT}/regional.grid.b     regional.grid.b
      /bin/ln -s ${RT}/depth_IASb0.50_03.a regional.depth.a
      /bin/ln -s ${RT}/depth_IASb0.50_03.b regional.depth.b
      touch   ${R}/archv.${ydh}.b
      /bin/rm ${R}/archv.${ydh}.[ab]
      ~/hycom/ALL/plot/src/trim_archv <<E-o-D
${R}/archv.${ydh}_L22.b
${R}/archv.${ydh}.b
000	'iexpt ' = experiment number x10 (000=from archive file)
  0	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
 93	'idm   ' = longitudinal array size
 65	'jdm   ' = latitudinal  array size
 22	'kdmold' = original number of layers
 21	'kdmnew' = target   number of layers
 25.0	'thbase' = reference density (sigma units)
 19.50	'sigma ' = layer  1  density (sigma units)
 20.25	'sigma ' = layer  2  density (sigma units)
 21.00	'sigma ' = layer  3  density (sigma units)
 21.75	'sigma ' = layer  4  density (sigma units)
 22.50	'sigma ' = layer  5  density (sigma units)
 23.25	'sigma ' = layer  6  density (sigma units)
 24.00	'sigma ' = layer  7  density (sigma units)
 24.70	'sigma ' = layer  8  density (sigma units)
 25.28	'sigma ' = layer  9  density (sigma units)
 25.77	'sigma ' = layer 10  density (sigma units)
 26.18	'sigma ' = layer 11  density (sigma units)
 26.52	'sigma ' = layer 12  density (sigma units)
 26.80	'sigma ' = layer 13  density (sigma units)
 27.03	'sigma ' = layer 14  density (sigma units)
 27.22	'sigma ' = layer 15  density (sigma units)
 27.38	'sigma ' = layer 16  density (sigma units)
 27.52	'sigma ' = layer 17  density (sigma units)
 27.64	'sigma ' = layer 18  density (sigma units)
 27.74	'sigma ' = layer 19  density (sigma units)
 27.82	'sigma ' = layer 20  density (sigma units)
 27.88	'sigma ' = layer 21  density (sigma units)
E-o-D
#
# --- reduce any resulting 22-layer archive files (*_L22.[ab]) to
# --- 21-layers using ALL/archive/src/trim_archv.  COMPLETED.
#
    endif
#
# --- add nested archive to tar_list.
#
    touch    ${R}/archv.${ydh}.b
    if (-z ${R}/archv.${ydh}.b) then
      echo "missing archive file: " ${R}/archv.${ydh}.b
      exit (2)
    else
      /bin/rm ${R}/archv.${ydh}_L22.[ab]
      echo archv.${ydh}.a >> tar_list.tmp
      echo archv.${ydh}.b >> tar_list.tmp
    endif
end
/bin/rm nest_list.tmp
#
# --- put nested archives in a tar bundle.
#
/bin/mv tar_list.tmp $R
cd $R
touch   archv_${Y}${A}.tar
/bin/rm archv_${Y}${A}.tar*
( cat tar_list.tmp | xargs tar cvf archv_${Y}${A}.tar ) >& archv_${Y}${A}.tar.lis
cat archv_${Y}${A}.tar.lis
#
# --- delete the individual nested archives.
#
cat tar_list.tmp | xargs /bin/rm
/bin/rm tar_list.tmp
