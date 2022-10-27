#!/usr/bin/gnuplot
set terminal pngcairo size 350,262 enhanced font 'Verdana,10'

filename = 'porazdelitev.txt'
unset key

system('mkdir -p proba')

stats filename name 'porazdelitev' nooutput

set xrange [0:960]
set yrange [0:500]

set xtics format ''
unset xtics
set ytics format ''
unset ytics

round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)

set style line 1 lc rgb '#8da0cb' lt 7 lw 2 ps .5 

#Make sure there's no blank line at the end to use STATS_blank !!
do for [i=1:int(porazdelitev_records)] {
    set output sprintf('proba/porazdelitev%03.0f.png',i)
    plot filename every ::0::i u 1:2 with points ls 1 notitle
}

# set output sprintf('images/porazdelitev1024.png')
# plot "porazdelitev1024.txt" u 1:2 with points ls 1 notitle


unset term
set term dumb
unset out
set out '/dev/null'
