#!/usr/bin/gnuplot
#set term cairolatex colour solid rounded noheader # size 12.5cm,10cm
#set out "graf3.tex"    #3 grafi za vsak plot posebej

#set key box top center width 1 height 0.1 Left reverse invert samplen 1
#set key opaque

set parametric
set isosamples 50,50
set hidden3d back

filename = 'fibonacci/932.txt'

#set grid ytics mytics x2tics mx2tics lc 'grey' dt 1
#set grid lc 'grey' dt 1

#set style line 1 lt 7 lw 2 ps .5
#set style fill transparent solid 0.5

set urange [-pi/2:pi/2]
set vrange [0:2*pi]

splot cos(u)*cos(v),cos(u)*sin(v),sin(u) w l lc 8 notitle,\
filename w points lt 7 lc 6 ps 1 notitle


#unset term
#set term dumb
#unset out
#set out '/dev/null'
