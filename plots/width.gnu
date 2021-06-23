reset
set terminal pdfcairo size 12, 6 font 'Libertinus Serif,20' lw 1
set output "growth-width.pdf"
set grid
set key box right bottom Left reverse height 0.5 width 0.5
set key opaque
LABEL = "y = 0.12868123x + 1.96174358"
set obj 1 rect at 25, 3 size char strlen(LABEL), char 2
set obj 1 front
set label at 25, 3 LABEL front center
set xlabel "time (ns)"
set ylabel "width (nm)"
plot\
    "../python/data/growth-width.dat" u 1:2 w l t "Simulated growth",\
    "../python/data/growth-linreg.dat" u 1:2 w l t "Linear fit"
    