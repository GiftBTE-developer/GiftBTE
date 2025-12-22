set xlabel 'X*'
set ylabel 'Temperature (K)'
plot "1e-8/TempLattice.dat" using ($2/1e-8):4 with p pt 4 ps 2 title "GiftBTE L=1e-8m", "1e-7/TempLattice.dat" using ($2/1e-7):4 with p pt 4 ps 2 title "GiftBTE L=1e-7m", "1e-6/TempLattice.dat" using ($2/1e-6):4 with p pt 4 ps 2 title "GiftBTE L=1e-6m", "analytic_201806/T1e-8.dat" using 1:($2+300) with l lt 1 lw 5 title "Analytical L=1e-8m", "analytic_201806/T1e-7.dat" using 1:($2+300) with l lt 2 lw 5 title "Analytical L=1e-7m", "analytic_201806/T1e-6.dat" using 1:($2+300) with l lt 3 lw 5 title "Analytical L=1e-6m"

