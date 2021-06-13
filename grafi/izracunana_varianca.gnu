#!/usr/bin/gnuplot
set term epslatex colour solid rounded noheader
set out "graf3.tex"

set key box top center width 1 height 0.1 Left reverse invert samplen 1
set key opaque

#set samples 10000
filename1 = '../thomson/varianca_povprecje_85.txt'
filename2 = '../fibonacci/varianca_povprecje_85.txt'
filename3 = '../random/varianca_povprecje_85.txt'

set style line 1 lt 7 lw 2 ps .5 lc rgb '#66c2a5'
set style line 2 lt 7 lw 2 ps .5 lc rgb '#fc8d62'
set style line 3 lt 7 lw 2 ps .5 lc rgb '#8da0cb'
#set style fill transparent solid 0.5

set title 'Izračunana varianca s povprečenjem za $N=85$'
set ylabel '$\sigma_N^2 (\vartheta)$'
set xlabel '$\vartheta/\pi$' 

plot filename1 w l ls 1 title 'Thomson',\
filename2 w l ls 2 title 'Fibonacci',\
filename3 w l ls 3 title 'Naključna'

unset term
set term dumb
unset out
set out '/dev/null'
