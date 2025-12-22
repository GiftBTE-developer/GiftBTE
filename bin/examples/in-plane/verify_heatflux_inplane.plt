set xlabel 'Y*'
set ylabel 'x-Heat Flux (W/m^2)'
plot "1e-8/HeatFlux.dat" using ($2/1e-8):4 with p pt 4 ps 2 title "GiftBTE L=1e-8m", "1e-7/HeatFlux.dat" using ($2/1e-7):4 with p pt 4 ps 2 title "GiftBTE L=1e-7m", "1e-6/HeatFlux.dat" using ($2/1e-6):4 with p pt 4 ps 2 title "GiftBTE L=1e-6m", "H1e-8.dat" using 1:2 with l lt 1 lw 5 title "Analytical L=1e-8m", "H1e-7.dat" using 1:2 with l lt 2 lw 5 title "Analytical L=1e-7m", "H1e-6.dat" using 1:2 with l lt 3 lw 5 title "Analytical L=1e-6m"

