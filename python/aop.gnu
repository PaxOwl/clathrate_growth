reset
set terminal pdfcairo size 12, 6 font 'Libertinus Serif,20' lw 1
set output "aop.pdf"
set grid
set key box right top Left reverse height 0.5 width 0.5
set key opaque
set xrange [0:10]
set xlabel "x (nm)"
set ylabel "aop"
plot\
    "rs_aop.dat" u 1:2 w l t "resampled aop",\
    