force/plot/README.force.plot:

Plots produced using the standard HYCOM 2-d horizontal field plot package.
All plots use NCAR graphics and are in logical space (i.e. every grid
cell is the same size on the plot).  

force/plot/README.force.plot	this file

force/plot/alias.src		aliases to simplify plot generation

force/plot/link.com		script to create softlinks
force/plot/regional.depth.a	softlink to bathymentry
force/plot/regional.depth.b	softlink to bathymentry
force/plot/regional.grid.a	softlink to grid
force/plot/regional.grid.b	softlink to grid
force/plot/src			softlink to latest source directory

force/plot/all.com		script to create all plots
force/plot/all.log		from csh all.com >& all.log
force/plot/coads_airtmp.IN	input for fieldproc
force/plot/coads_airtmp.ps	from fp2ps coads_airtmp.IN
force/plot/coads_precip.IN	input for fieldproc
force/plot/coads_precip.ps	from fp2ps coads_precip.IN
force/plot/coads_radflx.IN	input for fieldproc
force/plot/coads_radflx.ps	from fp2ps coads_radflx.IN
force/plot/coads_shwflx.IN	input for fieldproc
force/plot/coads_shwflx.ps	from fp2ps coads_shwflx.IN
force/plot/coads_vapmix.IN	input for fieldproc
force/plot/coads_vapmix.ps	from fp2ps coads_vapmix.IN
force/plot/coads_wndspd.IN	input for fieldproc
force/plot/coads_wndspd.ps	from fp2ps coads_wndspd.IN

The aliases in alias.src (source alias.src) can be used to simplify the
generation of plots.  For example,  fp2ps coads_airtmp.IN will generate
coads_airtmp.log and coads_airtmp.ps.

The plot program, fieldproc, is completely region independent and can
display the full domain or for a sub-region.  The number of plots per 
frame, spacing of latitude/longitude labels and grid lines, and the 
location of the contour label and color bar are all now specified at 
run-time.

The location of the sub-region is set at run time by specifying the 
location (iorign,jorign) on the full grid of (1,1) on the subregion 
grid and its size (idmp,jdmp).  For the full region, iorign=jorign=1 
and idmp=jdmp=0.  All other input parameters can be the same for the 
full region and a subregion.

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
