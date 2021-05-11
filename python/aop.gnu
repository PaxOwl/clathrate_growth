reset
set terminal pdfcairo size 12, 6 font 'Libertinus Serif,20' lw 1
set output "aop.pdf"
set grid
set key box right top Left reverse height 0.5 width 0.5
set key opaque
set xlabel "x (nm)"
set ylabel "aop"
plot\
    "aop_moy.dat" u 1:2 w l t "aop",\
    