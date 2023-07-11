set xlabel 'X*'
set ylabel 'T*'
plot "2e-7/TempLattice.dat" using ($2/2e-7):($4*8*150/1e15/2e-7/2e-7) with p pt 4 ps 2 title "GiftBTE L=2e-7m", "1e-6/TempLattice.dat" using ($2/1e-6):($4*8*150/1e15/1e-6/1e-6) with p pt 4 ps 2 title "GiftBTE L=1e-6m", "1e-5/TempLattice.dat" using ($2/1e-5):($4*8*150/1e15/1e-5/1e-5) with p pt 4 ps 2 title "GiftBTE L=1e-5m", "05.dat" using 1:2 with l lt 1 lw 5 title "Cao et al. L=2e-7m", "01.dat" using 1:2 with l lt 2 lw 5 title "Cao et al. L=1e-6m", "001.dat" using 1:2 with l lt 3 lw 5 title "Cao et al. L=1e-5m"

