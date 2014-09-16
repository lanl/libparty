set terminal svg size 800,600 enhanced fname "/usr/share/fonts/dejavu/DejaVuSansMono.ttf" fsize 18 dashed
set output "nnn.svg"
set style line 11 lc rgb '#000000' lt 1
set border 3 back ls 11
set tics nomirror
set style line 12 lc rgb '#808080' lt 3 lw 1
set grid back ls 12
set rmargin 5
set datafile separator ","
set logscale xy
set key bottom right
f1(x) = a1*x*log(x) +d1
d1=0.0001
fit f1(x) 'nnaos.txt' using 1:2 via a1, d1
f2(x) = a2*x*log(x) +d2
d2=0.0001
fit f2(x) 'nnsoa.txt' using 1:2 via a2, d2
f3(x) = a3*x*log(x) +d3
d3=0.0001
fit f3(x) 'nnaosoa.txt' using 1:2 via a3, d3
set style line 1 lt 1 lc rgb "#A00000" lw 3 pt 1
set style line 2 lt 1 lc rgb "#00A000" lw 3 pt 6
set style line 3 lt 1 lc rgb "#5060D0" lw 3 pt 2
set style line 4 lt 2 lc rgb "#A00000" lw 3 pt 1
set style line 5 lt 2 lc rgb "#00A000" lw 3 pt 6
set style line 6 lt 2 lc rgb "#5060D0" lw 3 pt 2

set xlabel "Number of Particles"
set ylabel "Seconds"
set yrange [.0001:]
set xrange [100:]
set mxtics 10
plot 'nnaos.txt' title "AOS" with lp ls 1, \
'nnsoa.txt' title "SOA" with lp ls 2, \
'nnaosoa.txt' title "AOSOA" with lp ls 3 , \
f1(x) title "AOS NlogN Fit" ls 4, \
f2(x) title "SOA NlogN Fit" ls 5, \
f3(x) title "AOSOA NlogN Fit" ls 6

