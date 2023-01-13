#!/usr/bin/gnuplot
set term epslatex colour solid rounded noheader
set size ratio -1
set out "porazdelitev_1024.tex"

set key box top left width -0.5 height 0.1 Left reverse invert samplen 0.5
set key opaque
unset key

#set samples 10000
filename = 'porazdelitev_1024.txt'

set style line 1 lt 7 lw 3 ps 1.5 lc rgb '#043565'
set style line 2 lt 7 lw 3 ps 1.5 lc rgb '#8f2d56'
set style line 3 lt 7 lw 3 ps 1.5 lc rgb '#218380'
set style line 4 lt 7 lw 3 ps 1.5 lc rgb '#ffc43d'
set style line 5 lt 7 lw 3 ps 1.5 lc rgb '#fb9f89'

unset xlabel
unset xtics 
unset ylabel
unset ytics 

plot filename u 1:2 w p ls 1 title 'tocke'


unset term
set term dumb
unset out
set out '/dev/null'