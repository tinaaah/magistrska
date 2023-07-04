#!/usr/bin/gnuplot
# set term cairolatex pdf color rounded linewidth 3 font 'cmp,12,m' size 14cm,10cm
set term cairolatex pdf color rounded linewidth 3 font 'cmp,12,m' size 13cm,8cm
set margins screen 0.1, screen 0.9, screen 0.1, screen 0.9
# set term tikz colour noheader size 12.5cm,10cm

#set samples 10000
filename1 = 'energy.txt'
filename2 = 'coord.txt'
filename3 = 'correlation.txt'
filename4 = 'rotation.txt'


set style line 1 lt 7 ps 1.5 lw 2 lc rgb "#f94144"
set style line 2 lt 7 ps 1.5 lw 2 lc rgb "#f3722c"
set style line 3 lt 7 ps 1.5 lw 2 lc rgb "#f8961e"
set style line 4 lt 7 ps 1.5 lw 2 lc rgb "#f9c74f"
set style line 5 lt 7 ps 1.5 lw 2 lc rgb "#90be6d"
set style line 6 lt 7 ps 1.5 lw 2 lc rgb "#43aa8b"
set style line 7 lt 7 ps 1.5 lw 2 lc rgb "#577590"

set border behind

# Energija ______________________________________________________ #
# set out "growth_a_01_energy.tex"
# set out "growth_a_01_stdenergy.tex"

# set xlabel '$a$'

# set ylabel '$u$'
# set ytics 20
  
# plot filename1 u 1:2 w l ls 6 notitle

# set ylabel '$\sigma_{u} / u$'
# set ytics 0.2

# plot filename1 u 1:3 w l ls 6 notitle 

# Stevilo sosedov _______________________________________________ #
# set out "growth_a_01_coord.tex"
# set out "growth_a_01_stdcoord.tex"

# set xlabel '$a$'

# set ylabel '$Z$'
# set ytics 1
  
# plot filename2 u 1:2 w l ls 6 notitle

# set ylabel '$\sigma_{Z} / Z$'
# set ytics 0.2
 
# plot filename2 u 1:3 w l ls 6 notitle 

# Orientacijska korelacijska funkcija ___________________________ #
# set out "growth_a_01_correlation.tex"
# set out "growth_a_01_stdcorrelation.tex"

# set xlabel '$r$'

# set ylabel '$\psi_2$'
# set ytics 0.2
 
# plot filename3 u 1:2 w l ls 6 notitle

# set ylabel '$\sigma_{\psi_2} / \psi_2$'
# set ytics 0.01
# 
# plot filename3 u 1:3 w l ls 6 notitle

# Sprejete in zavrnjene rotacije _________________________________ #
# set out "growth_a_01_rotations.tex"
set out "growth_a_01_stdrotations.tex"

set key box top right width 1 height 0.1 Left reverse invert samplen 0.5
set key opaque

set xlabel '$a$'
 
# set ylabel '$N$'
# set ytics 0.2
 
# plot filename4 u 1:2 w l ls 7 title "Sprejete",\
# filename4 u 1:4 w l ls 1 title "Zavrnjene"

set ylabel '$\sigma_N / N$'
set ytics 0.2
 
plot filename4 u 1:3 w l ls 7 title "Sprejete",\
filename4 u 1:5 w l ls 1 title "Zavrnjene"

unset out
