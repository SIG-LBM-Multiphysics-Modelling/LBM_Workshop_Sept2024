reset
set colorsequence classic
set terminal epslatex standalone color size 9cm,6cm
set out '2WedgePressureResultant.tex'
set yr[-0.:15]
set ytics 5
set xr[-0.:0.0041]
set xtics 0.001
set key top Left left samplen 2
set xlabel '$t\,[\mathrm{s}]$' offset 0,0.5
set ylabel '$P \, [\mathrm{kN/m}]$' offset -0.5,0
set pointsize 1.6
plot "data.txt" u 1:($3/1000) w l lc 0 lw 3 dt 1 title 'An.', "fsemi2_SPH.txt" u 1:($2/1000) w l lc 1 lw 3 dt 2 title 'SPH', "data.txt" u 1:($2/1000) w l lc 2 lw 3 dt 3 title 'FD', "Mom/data.txt" u 1:($2/1000) w l lc 3 lw 3 dt 4 title 'Mom' 
