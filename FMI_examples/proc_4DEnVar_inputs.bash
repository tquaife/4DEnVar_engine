#!/bin/bash -x

: '
This bash script reformats the data from FMI model output for input into the 4DEnVar code
'

#obs in years 5, 10 and 15
gawk '{print $5,$10,$15}' ensemblerun_1.dat > hx.dat.$$
#column 1 not needed
gawk '{print $2,$3,$4,$5,$6}' ensemble_22s_kg.dat > Xb.dat.$$
gawk -f flipit.awk Xb.dat.$$ > Xb.dat
#hx_bar
gawk '{for(i=1;i<=NF;i++)x[i]+=$i}END{for(i=1;i<=NF;i++)print x[i]/NR}'  hx.dat.$$ > hx_bar.dat
#transpose hx *after* calculating hx_bar
gawk -f flipit.awk hx.dat.$$ > hx.dat
#obs (y)
gawk '{print $2}' havainnot_t51015_1.dat > y.dat
#R, assuming diagonal. Note: need to convert to variances
gawk 'BEGIN{N=3}{for(i=1;i<=N;i++){if(i==NR){printf("%s ",$3*$3)}else{printf("0.00 ") }}print ""}' havainnot_t51015_1.dat > R.dat


#4DEnVar:
../4DEnVar Xb.dat hx.dat y.dat R.dat hx_bar.dat > OUT.$$

gawk '/===/{exit};1' OUT.$$ > xa.dat
gawk 'p;/===/{p=1}' OUT.$$ > Xa.dat

#plots:
python kdeplot.py

#tidy up
rm *.$$


