reset
set terminal pdfcairo size 12, 6 font 'Libertinus Serif,20' lw 1
set output "growth-cages.pdf"
set grid
set key box left center Left reverse height 0.5 width 0.5
set key opaque
set xlabel "Time (ns)"
set ylabel "Number of cages"
plot\
    "../python/data/growth-cages.dat" u 1:2 w l t "Small cages",\
    "../python/data/growth-cages.dat" u 1:3 w l t "Large cages",\
    "../python/data/growth-cages.dat" u 1:4 w l t "Interface cages",\
    "../python/data/growth-cages.dat" u 1:5 w l t "Irregular cages",\
    "../python/data/growth-cages.dat" u 1:6 w l t "Not in a cage"
    