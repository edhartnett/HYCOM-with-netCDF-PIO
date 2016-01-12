set terminal postscript default color solid
set output "sig.ps"
set title "Zonal Isopycnal Depths, ATLa2.00"
set xlabel "latitude"
set ylabel "depth (m)"
plot [-20:65] [0:60] \
    "sig_21.0" using 2:3 title "21.0" with lines, \
    "sig_21.5" using 2:3 title "21.5" with lines, \
    "sig_22.0" using 2:3 title "22.0" with lines, \
    "sig_22.5" using 2:3 title "22.5" with lines, \
    "sig_23.0" using 2:3 title "23.0" with lines, \
    "sig_23.5" using 2:3 title "23.5" with lines
plot [-20:65] [0:200] \
    "sig_23.5" using 2:3 title "23.5" with lines, \
    "sig_24.0" using 2:3 title "24.0" with lines, \
    "sig_24.5" using 2:3 title "24.5" with lines, \
    "sig_25.0" using 2:3 title "25.0" with lines, \
    "sig_25.5" using 2:3 title "25.5" with lines, \
    "sig_26.0" using 2:3 title "26.0" with lines
plot [-20:65] [0:1600] \
    "sig_26.0" using 2:3 title "26.0" with lines, \
    "sig_26.5" using 2:3 title "26.5" with lines, \
    "sig_27.0" using 2:3 title "27.0" with lines, \
    "sig_27.5" using 2:3 title "27.5" with lines

