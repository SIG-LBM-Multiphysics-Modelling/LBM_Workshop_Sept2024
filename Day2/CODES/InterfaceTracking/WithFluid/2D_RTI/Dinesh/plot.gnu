reset
set colorsequence classic
set terminal epslatex standalone color size 9cm,6cm
set out '2D_RTI_Re256.tex'
set xr[0:3]
set xtics 0.5
#set yr[-1.5:0]
set ytics 0.3
set key reverse bottom left Left samplen 1 width 1
set xlabel '$t/t_0$' offset 0,0.7
set ylabel '$ y^{\dagger}$' offset -0.,0
plot "Present.txt" u 1:($2/256) smooth bezier w l lc 0 dt 1 lw 3 title 'Present', "He.txt" u 1:($2+2) w l lc 1 dt 2 lw 3 title 'He $\textit{et al.}$',  "Dinesh.txt" u 1:($2+2) w l lc 2 dt 3 lw 3 title 'Dinesh $\textit{et al.}$'
