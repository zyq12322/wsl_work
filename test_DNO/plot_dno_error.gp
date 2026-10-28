set terminal pngcairo size 1200,700 enhanced font ",14"
set output "dno_error_plot.png"
set title "DNO (Airy wave): L2 error vs DNO FS order M"
set xlabel "M"
set ylabel "log_{10}( Relative L2 error )"
set grid
set key outside right top
set xrange [-0.5:15.5]

MbList = "0 1 2 3 4 5 6 7 8 9 10 11 12"
plot for [mb in MbList] "dno_error_results.txt" \
        using 1:($2 == mb + 0 ? log10($4) : 1/0) \
        with linespoints pt 7 lw 2 title "BT order Mb = ".mb

set output    # close the file so it is fully written
exit