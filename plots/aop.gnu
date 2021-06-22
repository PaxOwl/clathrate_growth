reset
set terminal pdfcairo size 12, 6 font 'Libertinus Serif,20' lw 1
set output "10melt-aop.pdf"
set grid
set key box right top Left reverse height 0.5 width 0.5
set key opaque
set xrange [0:10]
set xlabel "x (nm)"
set ylabel "AOP"
plot\
    "../python/data/10melt-300.dat" u 1:2 w l t "AOP",\
    