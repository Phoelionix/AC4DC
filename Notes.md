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

# Stability

Given the distinct features of the electron distribution, and the only obvious way to achieve stability being the computationally expensive grid points, a region-based approach with the grid points is necessary (at least, in lieu of numerical analysis of the basis splines). 

Divergences are more likely when an ionising process is occuring at a higher rate. Explicitly, the three "sharp" features of the distribution are the Maxwell-Boltzmann (MB) peak, the Auger peak, and the photoionisation peak. The MB peak's stability is apparently complicated by its peak energy going to 0 as t -> t_start. 



## Grid point region notes:
More detailed notes on each region follows. Energy ranges were intended to encapture Carbon, Nitrogen and Sulfur (and hence most proteins), but is based off Carbon. 6keV photons and a fluence of 2*10^7 J/cm^2 was used, but photon energy would only affect the photoionisation peak, and the fluence should not change the energy of the regions. 


eV      |
0-x     | empty, apart from point at 0. 
x       | CRITICAL, should be above  γ (current bounds are 1-5). The lower x is the more likely the simulation breaks/diverges there. 
    There are two regimes:
    x < γ:
        - simulation-breaking, cannot be remedied. Is either a separate effect or a regime where MB breaks the grid almost instantly
        Trying either using 3x the grid points or using ten times the time steps with x = 1 had a small but 
        impractical reduction on the divergence. 
    x>= γ:
        - small sad wiggle at x that is overtaken by MB relatively quick.  
        -Breakage more obviously due to MB peak being too sharp.
        - increasing time steps is helpful. For x = 5, going from 20k -> 100k steps moved divergence from t=-9.97 to t=-9.8552
        - Finess of grid points in 5-10 eV region doesn't seem to help, but it can make it break earlier with more of them.

10-50   | Sharp MB peak  - The sharpness of this peak is hidden by a log scale.
50-200  | Smooth for the early stages. However too few grid points here leads to crashes later. Worth running simulation in two parts with different grid points if this happens. (Still implementing) 
200-500 | Auger peak region, wider than MB peak on a linear scale, but sharper logarithmically.  (carbon peak 337 eV, I need to check why not ~ 270 eV)
500-3000/4500| transition region
3000-4000, 4500-6500| photoelectron delta-like peak, sharper than MB peak. Mainly made up of ionised core electrons for Carbon, but for heavier elements like Sulfur, the ionisation cross-section for higher shells is large enough that colder photoelectrons dominate. (To illustrate, in 1 fs for Carbon the median ionic charge is +1. For Sulfur, it's +8-+9.)


## Remedies for breaks in each region:
- Method considered a remedy if it increases the duration of stable simulation.
- Since spline overlap is sparse, it is probable that a divergence in one region is relatively independent of another.
- Of course, these regions are defined by myself, and there might be better classifications.  However, my goal 
  is to get what we need out of this simulation for templates, not to make it wiggle-proof. 

Region (ev) | Looks like | Remedy
x           | Divergence at x           | increase x!
x-10        | w/ higher grid point density than 10-50 region, very fine wiggles near end of region, but no wiggling (at first) below a point above x. Also had sudden disappearence with higher time steps.  | reduce number of points OR more time steps. (Possibility: Don't use a power law here?)
(Could be something else, combining this region with 10-50 leads to a wiggle minimum at 9.4 eV, presumably due to undershooting to reach up to MB)
leads to:
??          | small wiggle around 10 eV that keeps getting wider towards 0. | increase time steps.
10-50       | Wiggles in region that aren't divergent but make solution unstable    | increase grid fineness
50-200 | A continuous cycle of wiggles at MB peak that prevent MB from moving into this region (unless the intolerable stiff error is made low enough to catch it)    |  Simulate this period with finer grid points.
200-500     | Auger peak's width is unstable.  | increase grid fineness   (position was unstable a bit too when I varied some params, but haven't narrowed down cause)
4500-6500   | Lots of wiggles | Move to a linear scale where the entire peak can be seen. Only if it looks wiggly add grid points.


## TODO
Implement more sophisticated input settings that take in a csv/lists of a) the regions' energy boundaries, and b) a list of the power and number of grid points of each region    
