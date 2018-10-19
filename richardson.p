#gnuplot: graphical analysis of the Richardson method
reset
lb=-5.0
la=-1.0
N=10

T(x)=abs(x)>1.0 ? 0.5*(x+(x**2-1.0)**0.5)**N+0.5*(x-(x**2-1.0)**0.5)**N : cos(N*acos(x))

P(x)=T((2*x-la-lb)/(la-lb))/T((-la-lb)/(la-lb))

h(x)=2.0/(-lb-la+(lb-la)*cos((2.0*x-1)*pi/2.0/N))

set term qt 0
set xtics add (-2.0/4.0,-2.0/16.0)
set xrange [-7:0]
set ylabel "Amplitude"
set xlabel "Eigenvalue"
set arrow from -2.0/16.0,-0.2 to -2.0/16.0,1.0 nohead
set arrow from -2.0/4.0,-0.2 to -2.0/4.0,1.0 nohead
set arrow from -1.0,-0.2 to -1.0,1.0 nohead
plot P(x)

set term qt 1
set xrange [1:N]
unset arrow
set ylabel "Courant Factor"
set xlabel "Time Step Number"
plot h(x)
