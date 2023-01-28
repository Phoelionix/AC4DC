# Implementation Notes from earlier tests

## `tests/q_eii`
Plots `\sum_J Q_{xi; J, K}` and `\sum_{\eta} \Gamma_{\xi \to \eta; K}` against $K$, the energy-axis B-spline index
Good agreement is seen between the two (they should be the same) - however, this breaks if:
 - The grid is too sparse (n ~ 100 should be considered a minimum viable grid)
 - The minimum energy is 'too' high  (should be lower than smallest binding energy): This means that Dsigma in Qeii is not integrated over the full range of final states, causing disparities.

## `tests/test_abm.sh`
Automated test of the Adams integration routine at various orders. Test system is a harmonic oscillator.
 - Very high accuracies are shown to be achievable.

## `tests/rate_io_test`
Sanity check for rate IO functions.
Reads an existing EII.json, Auger.txt, Photo.txt, Fluor.txt and outputs to its second argument.
On success, it should satisfy 
` bin/tests/rate_io [infile] [outfile] && diff infile outfile` 

## `tests/basis_checker`
Populatrs a Distribution with a maxwellian of temperature T, then adds two deltas: a DeltaLike "implicit" spike at 3T and a DeltaSpike explicit source at 2T.

## `tests/spline_check`
Accepts same arguments as `BasisSet::set_parameters`, namely, `num_funcs, min, max, zero_degree_0, grid spacing`. Outputs all B-splines to stdout in a manner plottable by gnuplot

## `tests/sigma_test`
Verifies a few key properties of the sigma integrals.
1. Correctness: `∫dsigma(W|T) dW = 2 sigma(T)` (Integral may require large N to converge)
2. Symmetry: `dsigma(W|T) = dsigma(T-B-W | T)`

# Parameters
- Convergence is VERY sensitive to input parameters.
- Generally, wiggles in the 50-500eV region are attributable to QEII ad QTBR, but wiggles at lower energy are a sign of the influence of QEE.
- high-energy wiggles become less severe as the grid is made finer in that region.
- Low-energy wiggles are exacerbated by a finer grid - this is due to divergence of the derivative as x-> 0. Only affects QEE.

S.P. Notes:
# Parameters
- Computational time is proportional to the square of the total number of grid points.
