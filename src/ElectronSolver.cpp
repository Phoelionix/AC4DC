#include "ElectronSolver.h"
#include "HartreeFock.h"
#include "RateEquationSolver.h"
#include <fstream>

void PhotonFlux::set_parameters(double fluence, double fwhm){
    // The photon flux model
    // Gaussian model A e^-t^2/B
    B = fwhm*fwhm/(4*0.6931471806); // f^2/4ln(2)
    A = fluence/pow(Constant::Pi*B,0.5);
}
inline double PhotonFlux::operator()(double t){
    // Returns flux at time t
    return A*exp(-t*t/B);
};

Weight::Weight(size_t _size){
    size = _size;
    W = (double*) malloc(sizeof(double)*size*size); // allocate the memory
}

Weight::~Weight(){
    free(W);
}

double Weight::from(size_t idx){
    double x=0;
    for (size_t i = 0; i < size; i++) {
        x += operator()(i, idx);
    }
    return x;
}

double Weight::to(size_t idx){
    double x=0;
    for (size_t i = 0; i < size; i++) {
        x += operator()(idx, i);
    }
    return x;
}

ElectronSolver::ElectronSolver(const char* filename, ofstream& log) :
    MolInp(filename, log), pf(width, 100000*fluence)
{
    timespan = this->width*10;
}

void ElectronSolver::compute_cross_sections(std::ofstream& _log, bool recalc) {
    this->calc_rates(_log, recalc);
    hasRates = true;
    // Set up the rate equations!
    this->setup(state_type(Store, num_elec_points), this->timespan/num_time_steps);
    cout<<"Using timestep "<<this->dt<<" fs"<<endl;
};

void ElectronSolver::solve(){
    if (!hasRates) {
        std::cerr <<
"No rates have been loaded into the solver. \
Use ElectronSolver::compute_cross_sections(log)\n" << endl;
        return;
    }
    bool converged = false;
    // TODO: repeat this until convergence.
    this->iterate(-timespan/2, this->num_time_steps); // Inherited from ABM
}

//
//  !TODO
// void ElectronSolver::saveFree(const std::string&fname){
//     ofstream f;
//     f.open(fname);
//     f << "# Free electron dynamics"<<endl;
//     f << "# Energies \t Density [ UNITS ]" <<endl;
//     for (size_t i = 0; i < count; i++) {
//         f <<
//     }
//
// }

// The Big One: Incorporates all of the right hand side to the global
// d/dt P[i] = \sum_i=1^N W_ij - W_ji P[j]
// d/dt f = Q_B[f](t)
void ElectronSolver::sys(const state_type& s, state_type& sdot, const double t){
    for (size_t a = 0; a < s.atomP.size(); a++) {

        // create an appropriately-sized W[i][j]
        Weight W(s.atomP[a].size());

        double J = pf(t); // photon flux in [TODO: UNITS] units
        for ( auto& pht : Store[a].Photo) {
            W(pht.to, pht.from) = pht.val*J;
        }

        const state_type::bound_t& P = s.atomP[a];
        state_type::bound_t& Pdot = sdot.atomP[a];

        // Compute the changes in P
        for (size_t i = 0; i < P.size(); i++) {
            Pdot[i] = 0;
            for (size_t j = 0; j < i; j++) {
                Pdot[i] += W(i, j) * P[j] - W(j, i) * P[i];
            }
            // Avoid j=i
            for (size_t j = i+1; j < P.size(); j++) {
                Pdot[i] += W(i, j) * P[j] - W(j, i) * P[i];
            }
        }
    }
}
