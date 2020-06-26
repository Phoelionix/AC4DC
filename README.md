# AC4DC
Atomic/plasma code for XFEL damage simulations.

The code simulates the changes in electronic structure of a atoms in bulk of a micro-crystal or amorphous media as it is illuminated by an X-ray free-electron laser (XFEL) pulse. An independent atom approximation is adopted to represent a molecule as a collection of independent atoms, submerged into a bath of electrons that are formed as a result of ionization.

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

### Executable descriptions

`ac4dc`, the original version of the program, handles the electrons in the 'instant-thermalisation' approximation. Call with `ac4dc input_file(.txt|.inp)`, using `.inp` for atom data files and `.txt` for molecular specifications.

`rates`, the newer version, relaxes assumptions 4-6 and uses the approach of Leonov et al. to simulate the electron plasma dynamics.


### Output format

In the Xsections folder, the program stores a series of data files for Auger, Fluorescence and Photoionisation cross sections in the as `.txt`, in the format
`[ rate value ] [ index to ] [ index from ][ energy of transition ]`
Corresponding lines in the `index.txt` file explain which configuration the to/from indices correspond to.

EII parameters are stored in "sort-of-json" format - please note that the program does NOT use a robust json parser, and is sensitive to line formatting.


### TODO

1. Include some comments in the Xsection output file columns themselves
2. Electron plasma simulation part: sys
3. Move printing code to `state_type` namespace
4. Optimise with static arrays (if appropriate)

### Bibliography:

+ A. Kozlov and H. M. Quiney, _Comparison of Hartree-Fock and local density exchange approximations for calculation of radiation damage dynamics of light and heavy atoms in the field of x-ray free electron laser_, Phys. Scripta **94**, 075404 (2019). DOI: [10.1088/1402-4896/ab097c](https://doi.org/10.1088/1402-4896/ab097c)
+ C. P. Bhalla, N. O. Folland, and M. A. Hein, _Theoretical K-Shell Auger Rates, Transition Energies, and Fluorescence Yields for Multiply Ionized Neon_, Phys. Rev. A **8**, 649 (1973). DOI: [10.1103/PhysRevA.8.649](https://doi.org/10.1103/PhysRevA.8.649)
+ O. Yu. Gorobtsov, U. Lorenz, N. M. Kabachnik, and I. A. Vartanyants, _Theoretical study of electronic damage in single-particle imaging experiments at x-ray free-electron lasers for pulse durations from 0.1 to 10 fs_, Phys. Rev. E **91**, 062712 (2015). DOI: [10.1103/PhysRevE.91.062712](https://doi.org/10.1103/PhysRevE.91.062712)
+ W. R. Johnson, _Atomic Structure Theory: Lectures on Atomic Physics_, (Springer, Berlin, Heidelberg, 2007). DOI: [10.1007/978-3-540-68013-0](https://doi.org/10.1007/978-3-540-68013-0)
+ Leonov, A., D. Ksenzov, A. Benediktovitch, I. Feranchuk, and U. Pietsch, _Time Dependence of X-Ray Polarizability of a Crystal Induced by an Intense Femtosecond X-Ray Pulse._ IUCrJ **1**, no. 6 (2014): 402–17. DOI: [10.1107/S2052252514018156](https://doi.org/10.1107/S2052252514018156)
