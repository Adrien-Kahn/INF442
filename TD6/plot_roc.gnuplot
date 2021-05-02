unset label
set title "Courbe ROC"
set xlabel "Taux de fausses alarmes"
set ylabel "Taux de detection"
set pointsize 2
plot "knn_roc.out" with linespoints notitle linewidth 3
pause -1
