#
set echo
#
# --- test of curly vectors.
#
setenv NCARG_GKS_PS 013y021w_fsd+cur.ps
setenv OVERLAY VECTOR
#
./src/fieldtest <<'E-o-D'
../expt_01.3/data/archv.0021_016_00_fsd.a
../expt_01.3/data/archv.0021_016_00_3du.a
../expt_01.3/data/archv.0021_016_00_3dv.a
ATLb2.00
 57	'idm   ' = longitudinal array size
 52	'jdm   ' = latitudinal  array size
  1	'nperfr' = number of horizontal plots per frame
 10	'lalolb' = spacing of latitude/longitude labels
 10	'lalogr' = spacing of latitude/longitude grid over land (<0 land+sea)
  4	'loclab' = location of the contour label (1=upr,2=lowr,3=lowl,4=upl)
 11	'locbar' = location of the color bar     (1[0-4]=vert,2[0-4]=horiz)
  5	'kpalet' = palete (0=none,1=pastel,2=sst,3=gaudy,4=2tone,5=fc,6=ifc)
  1	'i_th  ' = draw a vector in only every i_th column
  1	'j_th  ' = draw a vector in only every j_th row
  1	'iorign' = i-origin of plotted subregion
  1	'jorign' = j-origin of plotted subregion
  0	'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
  0	'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
  1	'nrec  ' = next record to plot    (arbitrary order, <0 to end)
  1	'nrec2 ' = next record to overlay (arbitrary order)
FSD (cm) and V.1 - Jan Year 20
100.0	'qscale' = scale factor for plot
  2.0	'qq    ' = contour interval (<0 no plot; 0 from field)
  0.0	'center' = central contoured value
100.0	'vscale' = scale factor for vector (convert to cm/s)
 10.0	'vrefmx' = velocity plot maximum (streamline vectors)
  0.1	'vrefpl' = velocity plot maximum plot length (NDC: 0. to 1.)
  1	'nrec  ' = next record to plot    (arbitrary order, <0 to end)
  6	'nrec2 ' = next record to overlay (arbitrary order)
FSD (cm) and V.6 - Jan Year 20
100.0	'qscale' = scale factor for plot
  2.0	'qq    ' = contour interval (<0 no plot; 0 from field)
  0.0	'center' = central contoured value
100.0	'vscale' = scale factor for vector (convert to cm/s)
 10.0	'vrefmx' = velocity plot maximum (streamline vectors)
  0.1	'vrefpl' = velocity plot maximum plot length (NDC: 0. to 1.)
 -1	'nrec  ' = next record to plot (monotonicly increasing, <0 to end)
'E-o-D'
