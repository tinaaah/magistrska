#!/usr/bin/gnuplot
#set term epslatex colour solid rounded noheader
# graf od 25 do 28 luna | od 30 do 33 wolf | od 34 do 37 borza
#set out "graf.tex"

set key box top left width -0.5 height 0.1 Left reverse invert samplen 0.5
set key opaque

#set samples 10000
filename = 'centers.txt'

set style line 1 lt 7 lw 3 ps 1.5 lc rgb '#ef8a62'
set style line 2 lt 7 lw 3 ps 1.5 lc rgb '#67a9cf'
set style line 3 lt 7 lw 3 ps 1.5 lc rgb '#66c2a5'
set style line 4 lt 7 lw 3 ps 1.5 lc rgb '#998ec3'


plot filename u 1:2 w p ls 1 title 'tocke'


#unset term
#set term dumb
#unset out
#set out '/dev/null'