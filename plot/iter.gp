set terminal png enhanced font "arial,15" fontscale 1.0 size 1000, 1000 
set output "graphs/iter/iter_graphit.png"

set xlabel "Taille des données"
set ylabel "Temps en secondes"
set style data line
plot    "graphs/iter/iter_graphit.dat" using 1 title 'Richardson alpha',\
        "graphs/iter/iter_graphit.dat" using 2 title 'Jacobi',\
        "graphs/iter/iter_graphit.dat" using 3 title 'Gauss-Seidel'