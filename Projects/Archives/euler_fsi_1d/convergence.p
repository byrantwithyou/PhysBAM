plot "Standard_Tests/gnuplot_data.dat" using 1:2
f(x)=m*x+a
fit f(x) "Standard_Tests/gnuplot_data.dat" using 1:2 via m,a
replot f(x)
