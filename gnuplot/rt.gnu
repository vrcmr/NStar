#!/usr/bin/gnuplot 


#### comando para converter eps para pdf ps2pdf -dEPSFitPage file.eps file.pdf

reset
#set term postscript enhanced color
set terminal pngcairo crop
#set term png size 1149, 861 crop

set grid
#set xlabel "X"
#set ylabel "Y"
set xrange [-0.25:0.25]
set yrange [-0.75:0.75]
set zrange [*:*]
set cbrange [0.9:2.1]
unset clabel
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

nome = "../data/rt.0001.dat"

set view map
#set contour base
#set cntrparam levels 20
set palette rgbformulae 33,13,10
set cntrparam levels incremental 1.05,0.1,5


set output '../imagens/dens_rt_0001.png'
# set title "WENO5Z - RF"
splot nome u 1:2:3 w l lw 0.2 lc "black" notitle
#splot nome u 1:2:3 w l ls 7 palette notitle

#set cbrange [2.2:2.65]
#set output '../imagens/press_rt_0200.png'
#splot nome u ($1-0.25):($2-0.75):6 w l t "pressao"
