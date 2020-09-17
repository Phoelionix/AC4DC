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
Accepts same arguments as `BasisSet::set_parameters`, namely, `num_funcs, min, max, zero_degree, grid spacing`. Outputs all B-splines to stdout in a manner plottable by gnuplot

## `tests/sigma_test`
Verifies a few key properties of the sigma integrals.
1. Correctness: `∫dsigma(W|T) dW = 2 sigma(T)` (Integral may require large N to converge)
2. Symmetry: `dsigma(W|T) = dsigma(T-B-W | T)`

# What we Know
## Properties of QEII
- Satisfies `QEII[xi][J][K] =0` if `J-K > 3`, i.e. it looks like it ought to.
- Passes `check_qeii,` 
- Update - Now looks sensible. Cause of earlier problems: Unknown.

## Properties of QTBR
- Satisfies `\sum_J QTBR[xi][J][K][L] = \sum_eta \Gamma[xi->eta][K][L]` if J^2
   + K^2 > const, but not elsewhere.
- Symmetric in K and L (at least when summed over J)
- Low energy behavious is suspicions in Q: All K, L configurations should be
  able to three-body recombine.
- Hard edge present in Qtbr (as well as its complicated implementation) makes
  it more likely to be wrong.
- Hard edge depends on B and maximum energy. Corresponds to the line `ep + s = Emax - B/2.
- If ep + s < Emax - B/2, corresponds exactly to Gamma (as expected)
- Conclusion: Problem is TBR into regions with E > E_max. Sliced distribution
  represents the loss of electrons due to pairs being upconverted above the
energy ceiling.
  - Solution A: "she'll be right", leave the code as-is and rely on the
    vanishing of the distribution at large E to take care of it. Ensure that
max energy is at least twice as large as photoelectron peak.
	- Solution B: Fictional boundary. Enforces agreement between summed gamma and
	  Q rates by restricting the domain of the defining integrals to processes
    that are below the energy ceiling. Better matter conservation.

