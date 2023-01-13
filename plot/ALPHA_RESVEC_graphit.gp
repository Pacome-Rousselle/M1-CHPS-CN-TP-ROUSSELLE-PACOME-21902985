set terminal png enhanced font "arial,15" fontscale 1.0 size 1000, 1000 
set output "graphs/iter/ALPHA_RESVEC_graphit.png"

set xlabel "Norme du résidu"
set ylabel "Itération"
set style data line
plot  "graphs/iter/ALPHA_RESVEC_graphit.dat" notitle