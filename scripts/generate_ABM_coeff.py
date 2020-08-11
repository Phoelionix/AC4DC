#!/bin/python3
import numpy as np
from numpy.linalg import solve

# Uses the work of the paper
# https://www.researchgate.net/publication/249655976_A_Matrix_System_for_Computing_the_Coefficients_of_the_Adams_Bashforth-Moulton_Predictor-Corrector_formulae
# by Baba Seidu

MAX_ORDER = 15

f= open('src/Adams_arrays.h', 'w')
f.write('''
#ifndef ADAMS_ARRAY_CXX_H
#define ADAMS_ARRAY_CXX_H
// This file was written by scripts/generate_ABM_coeff.py gang
// Do not edit it manually, unless you really know what you're doing.

namespace AdamsArrays{
''')

f.write('''
const int MAX_ADAMS_ORDER = {0};
'''.format(MAX_ORDER))

f.write('const double AB_COEFF[%d][%d] = {' % (MAX_ORDER+1, MAX_ORDER+1))
for N in range(MAX_ORDER+1):
    if N == 0:
        f.write('\n\t{0}')
    elif N == 1:
        f.write(',\n\t{0}')
    else:
        # Adams-Bashforth predictor method
        A = np.array([[j**i for j in range(N)] for i in range(N)])
        C = np.array([(-1)**(i)/(i+1) for i in range(N)])
        B = solve(A,C)
        print('[Bashforth %d] obtained coeffs' % N, B)

        s = ',\n\t{'
        for i in range(len(B)):
            if i is not 0:
                s += ', '
            s += str(B[i])
        s += '}'
        f.write(s)

f.write('''
\n};

''')

f.write('const double AM_COEFF[%d][%d] = {' % (MAX_ORDER+1, MAX_ORDER+1))
for N in range(MAX_ORDER+1):
    if N == 0:
        f.write('\n\t{0}')
    elif N == 1:
        f.write(',\n\t{0}')
    else:
        # Adams-Moulton corrector method
        A = np.array([[(j-1)**i for j in range(N)] for i in range(N)])
        C = np.array([(-1)**(i)/(i+1) for i in range(N)])
        B = solve(A,C)
        print('[Moulton %d] obtained coeffs' % N, B)
        s = ',\n\t{'
        for i in range(len(B)):
            if i is not 0:
                s += ', '
            s += str(B[i])
        s += '}'
        f.write(s)

f.write('''
\n};

} // End of namespace\n

#endif
''')
f.close()
