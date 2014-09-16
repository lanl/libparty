set terminal svg size 800,600 enhanced fname "/usr/share/fonts/dejavu/DejaVuSansMono.ttf" fsize 18 dashed
set output "stepomp.svg"
set style line 11 lc rgb '#000000' lt 1
set border 3 back ls 11
set tics nomirror
set style line 12 lc rgb '#808080' lt 3 lw 1
set grid back ls 12
set rmargin 5
set datafile separator ","
#set logscale xy
set key top right
f1(x) = a1*(1/x*(1-d1) +d1)
fit f1(x) 'aosomp.txt' using 1:2 via a1, d1
f2(x) = a2*(1/x*(1-d2) +d2)
fit f2(x) 'soaomp.txt' using 1:2 via a2, d2
f3(x) = a3*(1/x*(1-d3) +d3)
fit f3(x) 'aosoaomp.txt' using 1:2 via a3, d3
set style line 1 lt 1 lc rgb "#A00000" lw 3 pt 1
set style line 2 lt 1 lc rgb "#00A000" lw 3 pt 6
set style line 3 lt 1 lc rgb "#5060D0" lw 3 pt 2
set style line 4 lt 2 lc rgb "#A00000" lw 3 pt 1
set style line 5 lt 2 lc rgb "#00A000" lw 3 pt 6
set style line 6 lt 2 lc rgb "#5060D0" lw 3 pt 2

set xlabel "Number of OpenMP Threads"
set ylabel "Seconds"
#set yrange [1000:]
set mxtics 10
#fittitle = sprintf("%.5f/N + %.5f", a1, d1)
#fittitle2 = sprintf("%.5f/N + %.5f", a3, d3)

fittitle = "AOS Amdahl's Law Fit"
fittitle2 = "SOA Amdahl's Law Fit"

plot 'aosomp.txt' title "AOS" with lp ls 1, \
'soaomp.txt' title "SOA" with lp ls 2, \
'aosoaomp.txt' title "AOSOA" with lp ls 3, \
f1(x) title fittitle ls 4, \
f3(x) title fittitle2 ls 6 #, \
#f3(x) title "AOSOA Linear Fit" ls 6

