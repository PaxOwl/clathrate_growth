reset
set terminal pdfcairo size 12, 6 font 'Libertinus Serif,20' lw 1
set output "100energy.pdf"
set grid
set key box left top Left reverse height 0.5 width 0.5
set key opaque
set datafile commentschars "#@&"
set xlabel "Time (fs)"
set ylabel "Potential Energy (kJ/mol)"
plot\
    "100energy270.xvg" u 1:2 w l t "270 K",\
    "100energy280.xvg" u 1:2 w l t "280 K",\
    "100energy290.xvg" u 1:2 w l t "290 K",\
    "100energy300.xvg" u 1:2 w l t "300 K",\
    "100energy310.xvg" u 1:2 w l t "310 K"