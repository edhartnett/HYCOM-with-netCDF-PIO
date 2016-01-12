plot/README.plot:

All plots use NCAR graphics and are in logical space (i.e. every grid
cell is the same size on the plot).  See ../archive for extracting
2-D and 3-D diagnostic fields into data files (several formats) suitable 
for plotting by other graphics packages.

plot/README.plot	this file
plot/README.plot.micom	how to plot MICOM archive files
plot/README.force.plot	how to plot scalar fields from any HYCOM ".a" file

plot/alias.src		aliases to simplify plot generation

plot/link.com		script to create softlinks
plot/regional.depth.a	softlink to bathymetry
plot/regional.depth.b	softlink to bathymetry
plot/regional.mask.b	softlink to land/sea mask (if any)
plot/regional.grid.a	softlink to grid
plot/regional.grid.b	softlink to grid
plot/src/		softlink to latest source directory

plot/010srf_020s.IN	text input for hycomproc (or hp), summer surface fields
plot/010srf_020s.ps	from  hp2ps 010srf_020s.IN
plot/010srf_mn.IN	text input for hycomproc (or hp), mean surface fields
plot/010srf_mn.ps	from  hp2ps 010srf_mn.IN
plot/010srf_sd.IN	text input for hycomproc (or hp), s.dev. surface fields
plot/010srf_sd.ps	from  hp2ps 010srf_sd.IN
plot/010y020s.IN	text input for hycomproc (or hp), summer
plot/010y020s.ps	from hp2ps 010y020s.IN
plot/010y021w.IN	text input for hycomproc (or hp), winter
plot/010y021w.ps	from hp2ps 010y020w.IN
plot/010y021w_sub.IN	text input for hycomproc (or hp), winter subregion
plot/010y021w_sub.ps	from hp2ps 010y021w_sub.IN

plot/depth.IN		text input for fieldproc (or fp), bathymetry
plot/depth.ps		from fp2ps depth.IN

The aliases in alias.src (source alias.src) can be used to simplify the
generation of plots.  For example,  hp2ps 010y020s.IN  will generate
010y020s.log and 010y020s.ps.

In this case, the only differences between summer and winter input file are:

diff 010y020s.IN 010y021w.IN
1c1
< ../expt_01.0/data/archv.0020_196_00.b
---
> ../expt_01.0/data/archv.0021_016_00.b

The first line identifies the archv file to plot and so must be different 
in every case.  The same plot program will also work with HYCOM 1.0
archive files, which are signaled by their filename (without an ".a" or
".b").

The input file allows very fine control over exactly what is plotted.
Not all layers need be included in the list of layer by layer plots.  
In fact, the same layer can appear more than once (giving fine control
over the plot order).  A negative number for 'kf' indicates that the 
layer number should be displayed on the plot, otherwise the layers
nominal isopycnal density is displayed.  If the color palete is 
multi-color (kpalet>1) and a positive contour interval is specified,
then the next input value is 'center' to identify the central value
of the color bar range.  The actual range then depends on the number 
of distinct colors in the palette (either 64 or 100):

      ipalet = 0      --  contour lines only, no color
      ipalet = 1      --  alternate pastel shading of contour intervals
      ipalet = 2      --  use canonical sst color palette  ( 64 intervals)
      ipalet = 3      --  use rainer's gaudy color palette (100 intervals)
      ipalet = 4      --  two-tone shades                  ( 64 intervals)
      ipalet = 5      --  NRL's 100 false color palette    (100 intervals)
      ipalet = 6      --  NRL's 100 inverted fc palette    (100 intervals)

Note that 'noisec' and 'noisec' (the number of z-sections to plot) must be 
followed by exactly the specified number of 'isec' and 'jsec' lines 
respectively.

The plot program, hycomproc, is completely region independent and can
display the full domain or for a sub-region.  It can also plot a "surface 
field only" archive file (by setting kk=1), and surface mean or
standard deviation files, see: 010srf*.IN and ../meanstd/README.meanstd.
The number of plots per frame, spacing of latitude/longitude labels and
grid lines, and the location of the contour label and color bar are all 
now specified at run-time.

The location of the sub-region is set at run time by specifying the location 
(iorign,jorign) on the full grid of (1,1) on the subregion grid and its
size (idmp,jdmp).  For the fill region, iorign=jorign=1 and idmp=jdmp=0.
All other input parameters can be the same for the full region and a 
subregion, but note that isec and jsec are w.r.t. the subregion rather 
than the full region.  See 010y021w_sub.IN.

The fieldproc program will plot scalar fields from an HYCOM 2.0 "*.a" file.
See depth.IN and README.force.plot.

Hycomproc can plot x-y surface fields, x-y layers, x-z slices and 
y-z slices, but fixed z-depth x-y plots require first archv2data3z 
(see ../archive) and then fieldproc.
