set terminal png enhanced font "arial,15" fontscale 1.0 size 1000, 1000 
set output "graphs/iter/errav_graphit.png"

set xlabel "Taille des données"
set ylabel "Erreur avant"
set style data line
plot    "graphs/iter/errav_graphit.dat" using 1 title 'Richardson alpha',\
        "graphs/iter/errav_graphit.dat" using 2 title 'Jacobi',\
        "graphs/iter/errav_graphit.dat" using 3 title 'Gauss-Seidel'