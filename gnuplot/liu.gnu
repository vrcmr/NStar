#!/usr/bin/gnuplot 


#### comando para converter eps para pdf ps2pdf -dEPSFitPage file.eps file.pdf

reset
#set term postscript enhanced color
set terminal pngcairo crop
#set term png size 1149, 861 crop

set grid
#set xlabel "X"
#set ylabel "Y"
set xrange [-.5:.5]
set yrange [-.5:.5]
set zrange [*:*]
set cbrange [0:5]
unset clabel
#unset key
set key outside
set size square

# set view 130, 10, 1, 1
#set samples 250, 250
set isosamples 20, 20
unset surface
set pm3d implicit at s
set pm3d scansforward


# cond. iniciais

nome = "../data/liu.0010.dat"

# set output '../imagens/dens3dliu.png'
# splot nome u ($1-.5):($2-.5):3 w l t "densidade"


set view map
set contour base
#set cntrparam levels 20
set palette rgbformulae 33,13,10
set cntrparam levels incremental 1.05,0.1,5


set output '../imagens/dens_liu_g.png'
# set title "WENO5Z - RF"
splot nome u ($1-.5):($2-.5):3 w l lw 0.2 lc "black" notitle
#splot nome u ($1-.5):($2-.5):3 w l ls 7 palette notitle

set output '../imagens/press_liu_g.png'
splot nome u ($1-.5):($2-.5):6 w l t "pressao"
