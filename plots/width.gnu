reset
set terminal pdfcairo size 12, 6 font 'Libertinus Serif,20' lw 1
set output "60width.pdf"
set grid
set key box right top Left reverse height 0.5 width 0.5
set key opaque
set xrange [0:25]
set xlabel "time (ns)"
set ylabel "width (Å)"
plot\
    "../python/data/60melt270-width.dat" u 1:2 w l t "270 K",\
    "../python/data/60melt280-width.dat" u 1:2 w l t "280 K",\
    "../python/data/60melt290-width.dat" u 1:2 w l t "290 K",\
    "../python/data/60melt300-width.dat" u 1:2 w l t "300 K"
    