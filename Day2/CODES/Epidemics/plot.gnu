reset
set colorsequence classic
set terminal epslatex standalone color size 9cm,6cm
set out 'SEIRD_Validation.tex'
set yr[-0.02:1.02]
set ytics 0.2
set xr[-2:100]
set xtics 25
set key center Left left samplen 2.5
set xlabel 'Days' offset 0,0.5
set ylabel 'Compartment fraction' offset 1,0
set pointsize 1.8
plot "dataBGK.txt" u 1:3 w l lc 0 lw 2 dt 1 title '$\hat{S}$', "dataFD.txt" u 1:3 every 500 w p lc 0 lw 2 dt 1 pt 5 notitle, "dataBGK.txt" u 1:4 w l lc 1 lw 2 dt 2 title '$\hat{E}$', "dataFD.txt" u 1:4 every 500 w p lc 1 lw 2 dt 2 pt 7 notitle,  "dataBGK.txt" u 1:5 w l lc 3 lw 2 dt 3 title '$\hat{I}$', "dataFD.txt" u 1:5 every 500 w p lc 3 lw 2 dt 2 pt 9 notitle, "dataBGK.txt" u 1:6 w l lc 4 lw 2 dt 4 title '$\hat{R}$', "dataFD.txt" u 1:6 every 500 w p lc 4 lw 2 dt 5 pt 11 notitle, "dataBGK.txt" u 1:7 w l lc 2 lw 2 dt 5 title '$\hat{D}$', "dataFD.txt" u 1:7 every 500 w p lc 2 lw 2 dt 8 pt 13 notitle

