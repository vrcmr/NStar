#!/usr/bin/gnuplot 


#### comando para converter eps para pdf ps2pdf -dEPSFitPage file.eps file.pdf

reset
#set term postscript enhanced color
set terminal pngcairo crop
#set term png size 1149, 861 crop

set grid
#set xlabel "X"
#set ylabel "Y"
set xrange [0.:7]
set yrange [0.:16]
set zrange [*:*]
#set cbrange [0.9:2.1]
#unset clabel
#unset key
set key outside
set size ratio 2



# set view 130, 10, 1, 1
#set samples 250, 250
set isosamples 20, 20
unset surface
set pm3d implicit at s
set pm3d scansforward


# cond. iniciais

nome = "../data/NS.0100.dat"

set view map
#set contour base
#set cntrparam levels 20
set palette rgbformulae 33,13,10
set cntrparam levels incremental 1.05,0.1,5

set yrange[5.5:7.5]
#set xrange[3.:5.]
#set cbrange [0:30]
set output '../imagens/press_NS_0100_total.png'
# set title "WENO5Z - RF"
splot nome u 1:2:6 w l lw 0.2 lc "black" notitle
#splot nome u 1:2:3 w l ls 7 palette notitle

set cbrange [0:700]
set output '../imagens/dens_NS_0060.png'
splot nome u 1:2:3 w l lw 0.2 lc "black" notitle
