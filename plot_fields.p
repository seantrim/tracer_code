#gnuplot
reset
nframes=1
aspect=1.0
comp=1
psize=1.0

unset key
set palette
set colorbox
set cbrange [*:*]
set xrange [0:aspect]
set yrange [0:1]
set size ratio 1.0/aspect
set pm3d map
unset xtics
unset ytics

set term qt 0
do for [i=0:nframes] {
n=10000+i
fname='T'.n
if (comp==1) {
set multiplot layout 2,3 rowsfirst
plot fname i 0 u 1:2:($3) w p pt 5 ps psize palette
plot fname i 0 u 1:2:($10) w p pt 5 ps psize palette
plot fname i 0 u 1:2:($4) w p pt 5 ps psize palette
plot fname i 0 u 1:2:($5) w p pt 5 ps psize palette
plot fname i 0 u 1:2:($6) w p pt 5 ps psize palette
plot fname i 0 u 1:2:($7) w p pt 5 ps psize palette}
else {
set multiplot layout 2,3 rowsfirst
plot fname i 0 u 1:2:($3) w p pt 5 ps psize palette
plot fname i 0 u 1:2:($4) w p pt 5 ps psize palette
plot fname i 0 u 1:2:($5) w p pt 5 ps psize palette
plot fname i 0 u 1:2:($6) w p pt 5 ps psize palette
plot fname i 0 u 1:2:($7) w p pt 5 ps psize palette
plot fname i 0 u 1:2:($10) w p pt 5 ps psize palette}
pause -1
unset multiplot}
