set datafile separator ","
set logscale xy
f1(x) = a1*x*log(x) +d1
fit f1(x) 'nntest' using 1:2 via a1, d1
f2(x) = a2*x*log(x) +d2
fit f2(x) 'nnomp' using 1:2 via a2, d2
f3(x) = a3*x*log(x) +d3
d3=10000
fit f3(x) 'nnomp2' using 1:2 via a3, d3
f4(x) = a4*x*log(x) +d4
fit f4(x) 'nnomp3' using 1:2 via a4, d4
f5(x) = a5*x*log(x) +d5
fit f5(x) 'nnomp4' using 1:2 via a5, d5
f6(x) = a6*x*log(x) +d6
fit f6(x) 'ompnomic' using 1:2 via a6, d6
f7(x) = a7*x*log(x) +d7
fit f7(x) 'ompnomic2' using 1:2 via a7, d7
plot 'nntest' with p, f1(x), 'nnomp' with p, f2(x), 'nnomp2' with p, f3(x), 'nnomp3' with p, f4(x), 'nnomp4' with p, f5(x), 'ompnomic' with p, f6(x), 'ompnomic2' with p, f7(x)


