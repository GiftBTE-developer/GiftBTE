set xlabel 'X*'
set ylabel 'T*'
plot "20nm/TempLattice.dat" using ($2/20e-9):4 with p pt 4 ps 2 title "GiftBTE L=20nm", "50nm/TempLattice.dat" using ($2/50e-9):4 with p pt 4 ps 2 title "GiftBTE L=50nm", "220nm/TempLattice.dat" using ($2/220e-9):4 with p pt 4 ps 2 title "GiftBTE L=220nm", "Tempcell_ref_20.dat" using ($1/20e-9):4 with l lt 1 lw 5 title "Ran et al. L=20m", "Tempcell_ref_50.dat" using ($1/50e-9):4 with l lt 2 lw 5 title "Ran et al. L=50m", "Tempcell_ref_220.dat" using ($1/220e-9):4 with l lt 3 lw 5 title "Ran et al. L=220m"

