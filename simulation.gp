set title "EnzymePFK strange Attractor "
set xlabel "X"
set ylabel "Y %"
set zlabel "Z %"
splot "simulation.dat" using 1:2:3 title 'Van Der Pol Attractor' with lines
