!/usr/bin/gnuplot -p

plot "data1.dat" using 1:2 title "b_0 = 70.0 (km/sec/Mpc), Lambda = 0.0 (1/m^2)" w l, "data2.dat" using 1:2 title "b_0 = 70.0 (km/sec/Mpc), Lambda = 1.25e-52 (1/m^2)" w l, "data3.dat" using 1:2 title "b_0 = 100.0 (km/sec/Mpc), Lambda = 0.0 (1/m^2)" w l, "data4.dat" using 1:2 title "b_0 = 100.0 (km/sec/Mpc), lambda = 1.25e-52 (1/m^2)" w l
set xlabel "Time (Billion Years Ago)"
set ylabel "Scale factor"
#set yrange [0:1.2]
#set xrange [0:1]
set title "phys3071 as06 melsom 42593249"
set term postscript color
set output "as06-problem2-melsom-42593249.ps"
replot
