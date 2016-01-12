relax/plot/README.relax.plot:

Plots produced using the standard HYCOM archive plot package from the
dummy archive version of each monthly climatology.  All plots use NCAR 
graphics and are in logical space (i.e. every grid cell is the same 
size on the plot).

relax/plot/README.relax.plot	this file

relax/plot/alias.src		aliases to simplify plot generation

relax/plot/link.com		script to create softlinks
relax/plot/regional.depth.a	softlink to bathymentry
relax/plot/regional.depth.b	softlink to bathymentry
relax/plot/regional.grid.a	softlink to grid
relax/plot/regional.grid.b	softlink to grid
relax/plot/src			softlink to latest source directory

relax/plot/all_010.com		script to create all plots for expt 1.0
relax/plot/all_010.log		from csh all_010.com >& all_010.log
relax/plot/010_jan_cs1.IN	input for hp, january  500m cross-sections
relax/plot/010_jan_cs1.ps	from hp2ps 010_jan_cs1.IN
relax/plot/010_jan_cs2.IN	input for hp, january 5500m cross-sections
relax/plot/010_jan_cs2.ps	from hp2ps 010_jan_cs2.IN
relax/plot/010_jan_hor.IN	input for hp, january horizontal plots
relax/plot/010_jan_hor.ps	from hp2ps 010_jan_hor.IN
relax/plot/010_jul_cs1.IN	input for hp, july     500m cross-sections
relax/plot/010_jul_cs1.ps	from hp2ps 010_jul_cs1.IN
relax/plot/010_jul_cs2.IN	input for hp, july    5500m cross-sections
relax/plot/010_jul_cs2.ps	from hp2ps 010_jul_cs1.IN
relax/plot/010_jul_hor.IN	nput for hp, july    horizontal plots
relax/plot/010_jul_hor.ps	from hp2ps 010_jul_hor.IN

The first line of the input text file identifies the dummy archive file
file to plot, for example:

ajax 78> diff 010_jan_hor.IN 010_jul_hor.IN
1c1
< ../010/relax.0000_016_00.a
---
> ../010/relax.0000_196_00.a

The aliases in alias.src (source alias.src) can be used to simplify the
generation of plots.  For example,  hp2ps 010_jul_hor.IN  will generate
010_jul_hor.log and 010_jul_hor.ps.

The input file allows very fine control over exactly what is plotted.
Not all layers need be included in the list of layer by layer plots.  
In fact, the same layer can appear more than once (giving fine control
over the plot order).  A negative number for 'kf' indicates that the 
layer number should be displayed on the plot, otherwise the layers
nominal isopycnal density is displayed.  Note that 'noisec' and 'noisec' 
(the number of z-sections to plot) must be followed by exactly the 
specified number of 'isec' and 'jsec' lines respectively.

If the color palete is multi-color (kpalet>1) and a positive contour 
interval is specified, then the next input value is 'center' to identify 
the central value of the color bar range.  The actual range then depends 
on the number of distinct colors in the palette (either 64 or 100):

      ipalet = 0      --  contour lines only, no color
      ipalet = 1      --  alternate pastel shading of contour intervals
      ipalet = 2      --  use canonical sst color palette  ( 64 intervals)
      ipalet = 3      --  use rainer's gaudy color palette (100 intervals)
      ipalet = 4      --  two-tone shades                  ( 64 intervals)
      ipalet = 5      --  NRL's 100 false color palette    (100 intervals)
      ipalet = 6      --  NRL's 100 inverted fc palette    (100 intervals)

The plot program, hycomproc, is completely region independent and can
display the full domain or a sub-region.  The number of plots per 
frame, spacing of latitude/longitude labels and grid lines, and the 
location of the contour label and color bar are all now specified at 
run-time.

The location of the sub-region is set at run time by specifying the location 
(iorign,jorign) on the full grid of (1,1) on the subregion grid and its
size (idmp,jdmp).  For the full region, iorign=jorign=1 and idmp=jdmp=0.
All other input parameters can be the same for the full region and a 
subregion, but note that isec and jsec are w.r.t. the subregion rather 
than the full region.
