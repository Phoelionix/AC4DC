# AC4DC
Atomic/plasma code for XFEL damage simulations.

The code simulates the changes in electronic structure of a atoms in bulk of a micro-crystal or amorphous media as it is illuminated by an X-ray free-electron laser (XFEL) pulse. An independent atom approximation is adopted to represent a molecule as a collection of independent atoms, submerged into a bath of electrons that are formed as a result of ionization.

### Versions

+ Version 1: Heuristic Values
  + Computes cross-sections for atomic processes: photoionisation, fluorescence, Auger decay and electron-impact ionisation
  + Solves dynamical rate equations for atomic state population evolution in time, assuming that electrons in the plasma follow a thermal distribution
+ Version 2: Totally Not Thermal
  + Allows an arbitrary form of the free-electron distribution. _f(E)_ interpolated over a basis of B-splines.
  + Concurrently solves for the time evolution of the free electrons and bound states, via the Boltzmann transport equation.

### Dependencies:

+ Eigen 3        (vector algebra library, http://eigen.tuxfamily.org/index.php?title=Main_Page)
+ Wigner Symbols (3-j and 6-j symbols, https://github.com/joeydumont/wignerSymbols)
+ OpenMP         (multi-threading, https://www.openmp.org/)

### Installation

Compatibility is only promised for Linux variants and macOS. (lack of Windows support is mainly due to the UNIX-style path assumptions used throughout. In principle, this may be corrected for by rewriting with `boost::filesystem`.)
Tested environments are:
+ Debian 10 (buster)
+ macOS 10.14.6 (Mojave)
This project is built with `make`, but does not have a `configure` script. Some fiddling may be required depending on your OS.
**Note for macOS users**  - Apple's standard `gcc` distribution does not support openmp, used for parallelisation. Code was compiled with `g++-9` from Homebrew. In principle, Homebrew llvm's `clang++` should also work, but has not been tested.

### Key Approximations in ac4dc:

1. All atoms interact with X-rays independent of each other.
2. At any time, an atom is described by an average-over-configuration state.
3. Time evolution of atoms is described by a system of coupled rate equation (perturbation theory).
4. Electron plasma is split into two components - energetic photo-electrons and cold secondary electrons.
5. Photo-electrons have delta-function distribution around time-dependent average kinetic energy. No three-body recombination (TBR) for photo-electrons, which can escape the sample.
6. Secondary electrons have Boltzmann distribution with time-dependent number density and temperature. They are assumed to be trapped at all times. Electron impact ionisation (EII) and TBR are determined by average number density of secondary electrons over the volume of a sample.

### Key Approximations in ac4dc2:

1. All atoms interact with X-rays independent of each other.
2. At any time, an atom is described by an average-over-configuration state.
3. Time evolution of atoms is described by a system of coupled rate equation (perturbation theory).
4. There are enough free electrons for a Boltzmann hydrodynamic treatment to be appropriate.
5. Electron plasma is isotropic and homogeneous in position space; isotropic in momentum space.

### Running AC4DC

Run with
`./ac4dc2 input/Neon.mol`

Flags: (after the inout file name)
- `-h`: Displays a brief help message
- `-s`: Authorises the program to (s)earch the `output` directory for previously computed atomic data. Note that there is currently no check performed to guarantee that these rates were calculated for the correct parameters, so this flag is not recommended outside of debugging. Otherwise, new parameters are calculated from scratch.

### Output format

In the Xsections folder, the program stores a series of data files for Auger, Fluorescence and Photoionisation cross sections in the as `.txt`, in the format
`[ rate value ] [ index to ] [ index from ][ energy of transition ]`
Corresponding lines in the `index.txt` file explain which configuration the to/from indices correspond to.

EII parameters are stored in "sort-of-json" format - please note that the program does NOT use a robust json parser, and is sensitive to line formatting.


### TODO

1. Try to eliminate numerical artefacts when running with coarse grids at low energy
4. Refactor to get rid of compiler warnings wherever possible
5. Implement 'output version control' for atomic parameters in storage: avoid unnecessary recalculation, guarantee recalculation if new input parameters are incompatible
  - Add methods to `Input.cpp` to enable reading/writing salient parameters to file, e.g. `output/C/run_2021-04-11/input.txt`
  - Add linear search implementation to input logic
6. Cmake build system
7. Restructure parameter input and rate output files to use JSON format
3. Optimise with static arrays - promote state_type to a N_FREE-dimensioned template for faster reads. (remains to be seen if this is a bottleneck)

### Bibliography:

+ A. Kozlov and H. M. Quiney, _Comparison of Hartree-Fock and local density exchange approximations for calculation of radiation damage dynamics of light and heavy atoms in the field of x-ray free electron laser_, Phys. Scripta **94**, 075404 (2019). DOI: [10.1088/1402-4896/ab097c](https://doi.org/10.1088/1402-4896/ab097c)
+ C. P. Bhalla, N. O. Folland, and M. A. Hein, _Theoretical K-Shell Auger Rates, Transition Energies, and Fluorescence Yields for Multiply Ionized Neon_, Phys. Rev. A **8**, 649 (1973). DOI: [10.1103/PhysRevA.8.649](https://doi.org/10.1103/PhysRevA.8.649)
+ O. Yu. Gorobtsov, U. Lorenz, N. M. Kabachnik, and I. A. Vartanyants, _Theoretical study of electronic damage in single-particle imaging experiments at x-ray free-electron lasers for pulse durations from 0.1 to 10 fs_, Phys. Rev. E **91**, 062712 (2015). DOI: [10.1103/PhysRevE.91.062712](https://doi.org/10.1103/PhysRevE.91.062712)
+ W. R. Johnson, _Atomic Structure Theory: Lectures on Atomic Physics_, (Springer, Berlin, Heidelberg, 2007). DOI: [10.1007/978-3-540-68013-0](https://doi.org/10.1007/978-3-540-68013-0)
+ Leonov, A., D. Ksenzov, A. Benediktovitch, I. Feranchuk, and U. Pietsch, _Time Dependence of X-Ray Polarizability of a Crystal Induced by an Intense Femtosecond X-Ray Pulse._ IUCrJ **1**, no. 6 (2014): 402–17. DOI: [10.1107/S2052252514018156](https://doi.org/10.1107/S2052252514018156)
