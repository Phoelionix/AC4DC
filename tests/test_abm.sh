#!/bin/bash

echo "============================================================"
echo "  Testing ODE method: SHO "
echo "Compiling tests/abm_verif.cpp"
g++ -g -std=c++11 -O0 tests/abm_verif.cpp -Isrc -o bin/tests/abm_verif
bin/tests/abm_verif 300 1e-6 5 > tmp/abm5_h6.txt
gnuplot -e "prefix='abm5_h6'" tests/abm_test.gnuplot
bin/tests/abm_verif 300 1e-6 3 > tmp/abm3_h6.txt
gnuplot -e "prefix='abm3_h6'" tests/abm_test.gnuplot
bin/tests/abm_verif 300 1e-3 8 > tmp/abm8_h3.txt
gnuplot -e "prefix='abm8_h3'" tests/abm_test.gnuplot
bin/tests/abm_verif 300 1e-4 8 > tmp/abm8_h4.txt
gnuplot -e "prefix='abm8_h4'" tests/abm_test.gnuplot
bin/tests/abm_verif 300 1e-6 8 > tmp/abm8_h6.txt
gnuplot -e "prefix='abm8_h5'" tests/abm_test.gnuplot
echo "Output saved to testOutput";
echo "============================================================"
