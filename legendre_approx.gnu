#!/usr/bin/gnuplot
#set term cairolatex colour solid rounded noheader # size 12.5cm,10cm
#set out "graf2.tex"    #3 grafi za vsak plot posebej

#set key box top center width 1 height 0.1 Left reverse invert samplen 1
#set key opaque

#set samples 10000
filename1 = 'proba_napaka.txt'
#filename1 = 'varianca/thomson_85.txt'
#filename2 = 'napake_legendre/random_85.txt'

#set grid ytics mytics x2tics mx2tics lc 'grey' dt 1
#set grid lc 'grey' dt 1

#set style line 1 lt 7 lw 2 ps .5
#set style fill transparent solid 0.5

plot filename1 w l lc -1 lw 2 notitle#,\
#filename2 w l lc 7 lw 2 notitle

#filename every :::0::6 w lp ls 1 lc rgb '#5066c2a5' notitle 'diskretni nivoji (sosedi)',\
#1/0 w p lt 5 ps 1.25 lc rgb '#508da0cb' title 'diskretni nivoji (poljubni)',\
#filename every :::7::13 w lp ls 1 lc rgb '#508da0cb' notitle 'diskretni nivoji (poljubni)',\
#1/0 w p lt 5 ps 1.25 lc rgb '#50fc8d62' title 'zvezni nivoji',\
#filename every :::14::20 w lp ls 1 lc rgb '#50fc8d62' notitle 'zvezni nivoji'


#unset term
#set term dumb
#unset out
#set out '/dev/null'
