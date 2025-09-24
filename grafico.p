set title 'Série Temporal'                       # plot title
# set title 'Caso Base ({/Symbol b} = 3.5)'                       # plot title
set xlabel 'Tempo'                              # x-axis label
set ylabel 'Pessoas Infectadas (Normalizado)'                          # y-axis label
# set ylabel 'Individuos'                          # y-axis label


# labels
# set label "boiling point" at 10, 212

# key/legend
# set key top right
# set key box
set key box vertical width 2 height 1 maxcols 1 spacing 1
set key at graph 0.95, 0.95
# set format y "%.0s*10^%T"; set ytics add ('0' 0)
# set format y "%2.0t*10^{%L}"


set datafile separator ',' 

plot "sir_rk4.csv" using 1:2 with lines title "S" lw 2, \
     'sir_rk4.csv' using 1:3 with lines title "I" lw 2, \
     'sir_rk4.csv' using 1:4 with lines title "R" lw 2, \

set terminal pngcairo
set output 'sir_rk4.png'
replot
unset output