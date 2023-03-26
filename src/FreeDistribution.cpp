/**
 * @file FreeDistribution.cpp
 * @brief @copybrief FreeDistribution.h
 * @details Grid points are synonymous w/ knot points for the B-splines. Note that this code expands 
 * on the original AC4DC model by introducing these grid points, allowing for an approximation of a 
 * continuous electron density distribution along the energy basis. 
 */
/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================*/

#include "FreeDistribution.h"
#include "Dipole.h"
#include "Constant.h"
#include "SplineIntegral.h"
#include <Eigen/StdVector>

// #define NDEBUG
// to remove asserts

// Initialise static things

size_t Distribution::size=0;  // modified by set_distribution now - S.P.
size_t Distribution::CoulombLog_cutoff=0;
double Distribution::CoulombDens_min=0;
SplineIntegral basis;

// Psuedo-constructor thing (Maybe not anymore... -S.P.)
void Distribution::set_basis(size_t step, GridSpacing grid_style, Cutoffs param_cutoffs, FeatureRegimes regimes, GridBoundaries elec_grid_regions){
    // Defines a grid of num_funcs points (if manual, thsi is the total number of free-electron grid points specified in the .mol file.)
    // where num_funcs is the number of non-boundary (i.e. "usable") splines/knots.
    //basis_history.push_back(SplineIntegral());
    basis.set_parameters(grid_style, elec_grid_regions,regimes);
    //basis_history.push_back(make_pair(step,basis));
    Distribution::size=basis.num_funcs;
    Distribution::CoulombLog_cutoff = basis.i_from_e(param_cutoffs.transition_e);
    Distribution::CoulombDens_min = param_cutoffs.min_coulomb_density;
    // cout<<"[ Free ] Estimating lnLambda based on first ";
    // cout<<CoulombLog_cutoff<<" points, up to "<<grid_style.transition_e<<" Ha"<<endl;    
    cout<<"[ Free ] Estimating lnLambda for energies below ";
    cout<<param_cutoffs.transition_e<<" Ha"<<endl;
    cout<<"[ Free ] Neglecting electron-electron below density of n = "<<CoulombDens_min<<"au^-3"<<endl;
}

void Distribution::set_distribution(vector<double> new_knot, vector<double> new_spline_factors) {
    f = new_spline_factors; 
    // Remove boundary knots
    new_knot = get_trimmed_knots(new_knot);
    size = new_knot.size();
    basis = new_knot;
}

// Adds Q_eii to the parent Distribution
void Distribution::get_Q_eii (Eigen::VectorXd& v, size_t a, const bound_t& P, const int threads) const {
    assert(basis.has_Qeii());
    assert(P.size() == basis.Q_EII[a].size());
    assert((unsigned) v.size() == size);
    for (size_t xi=0; xi<P.size(); xi++) {
        // Loop over configurations that P refers to
        double v_copy [size] = {0};
        #pragma omp parallel for num_threads(threads) reduction(+ : v_copy)
        for (size_t J=0; J<size; J++) {
            for (size_t K=0; K<size; K++) {
                v_copy[J] += P[xi]*f[K]*basis.Q_EII[a][xi][J][K];
            }
        }
        v += Eigen::Map<Eigen::VectorXd>(v_copy,size);
    }
}

/**
 * @brief 
 * @details d/dt f = Q_B[f](t)  
 * 
 * @param v 
 * @param a 
 * @param P d/dt P[i] = \sum_i=1^N W_ij - W_ji ~~~~ P[j] = d/dt(average-atomic-state)
 */
// Puts the Q_TBR changes in the supplied vector v
void Distribution::get_Q_tbr (Eigen::VectorXd& v, size_t a, const bound_t& P, const int threads) const {
    assert(basis.has_Qtbr());
    assert(P.size() == basis.Q_TBR[a].size());
    double v_copy [size] = {0}; 
    #pragma omp parallel for num_threads(threads) reduction(+ : v_copy) collapse(2)
    for (size_t eta=0; eta<P.size(); eta++) {          
        // Loop over configurations that P refers to
        for (size_t J=0; J<size; J++) {                   // size = num grid points
            for (auto& q : basis.Q_TBR[a][eta][J]) {   // Thousands of iterations for each J - S.P.
                 v_copy[J] += q.val * P[eta] * f[q.K] * f[q.L];
            }
        }
    }
    v += Eigen::Map<Eigen::VectorXd>(v_copy,size);
}

// Puts the Q_EE changes in the supplied vector v
void Distribution::get_Q_ee(Eigen::VectorXd& v, const int threads) const {
    assert(basis.has_Qee());
    double CoulombLog = this->CoulombLogarithm();
    // double CoulombLog = 3.4;
    if (CoulombLog <= 0) return; // skip calculation if LnLambda vanishes

    if (isnan(CoulombLog) || CoulombLog > 11.5) CoulombLog = 11.5;
    // double LnLambdaD = 0.5*log(this->k_temperature()/4/Constant::Pi/this->density());
    // if (isnan(LnLambdaD)) LnLambdaD = 11;
    // cerr<<"LnDebLen = "<<LnLambdaD<<endl;
    // A guess. This should only happen when density is zero, so Debye length is infinity.
    // Guess the sample size is about 10^5 Bohr. This shouldn't ultimately matter much.   /// Attention - S.P.
    double v_copy [size] = {0}; 
    #pragma omp parallel for num_threads(threads) reduction(+ : v_copy)  collapse(2)       
    for (size_t J=0; J<size; J++) {
        for (size_t K=0; K<size; K++) {
            for (auto& q : basis.Q_EE[J][K]) {
                 v_copy[J] += q.val * f[K] * f[q.idx] * CoulombLog;
            }
        }
    }
    v += Eigen::Map<Eigen::VectorXd>(v_copy,size);
}

void Distribution::get_Jac_ee(Eigen::MatrixXd& M) const{   // Unused currently -S.P.
    // Returns Q^p_qjc^q + Q^p_jrc^r
    assert(basis.has_Qee());
    M = Eigen::MatrixXd::Zero(size,size);
    double CoulombLog=3;
    // double CoulombLog = CoulombLogarithm(size/3);
    if (isnan(CoulombLog) || CoulombLog <= 0) return;
    for (size_t P=0; P<size; P++) {
        for (size_t Q=0; Q<size; Q++) {
            for (auto& q : basis.Q_EE[P][Q]) {
                // q.idx is  L
                M(P,q.idx) += q.val * f[Q] * CoulombLog;
                M(P, Q) += q.val * f[q.idx] * CoulombLog;
            }
        }
    }
}

/*
// Sets f_n+1 based on using a Newton-Rhapson-Euler stiff solver
// Not currently used
void Distribution::from_backwards_Euler(double dt, const Distribution& prev_step_dist, double tolerance, unsigned maxiter){

    Eigen::MatrixXd M(size, size);
    Eigen::MatrixXd A(size, size);
    Eigen::VectorXd v(size);
    Eigen::VectorXd old_curr_step(size);
    
    unsigned idx = 0;
    const unsigned MAX_ITER = 200;

    bool converged = false;

    Eigen::Map<const Eigen::VectorXd > prev_step (prev_step_dist.f.data(), size);
    Eigen::Map<Eigen::VectorXd > curr_step(f.data(), size);

    while (!converged && idx < MAX_ITER){
        
        old_curr_step = curr_step;
        v = Eigen::VectorXd::Zero(size);
        this->get_Q_ee(v, size/3); // calculates v based on curr_step
        
        // Newton iterator
        get_Jac_ee(M);
        A = dt * basis.Sinv(M) - Eigen::MatrixXd::Identity(size,size);
        
        curr_step = -A.fullPivLu().solve(dt*v + prev_step - curr_step);
        curr_step += prev_step;
        
        // Picard iterator
        // curr_step = prev_step + dt*v;

        double delta = fabs((curr_step-old_curr_step).sum());

        if (delta/curr_step.sum() < tolerance) converged = true;
        idx++;
    }
    if (idx == MAX_ITER) std::cerr<<"Max stiff-solver iterations exceeded"<<std::endl;
}
*/

// // Taken verbatim from Rockwood as quoted by Morgan and Penetrante in ELENDIF
// void Distribution::add_Q_ee(const Distribution& d, double kT) {
//     double density=0;
//     double CoulombLog = log(kT/(4*Constant::Pi*density));
//     double alpha = 2*Constant::Pi*sqrt(2)/3*CoulombLog;


// }

// ============================
// utility

// See add_density_destribution comments
// Expects T to have units of Ha
void Distribution::add_maxwellian(double T, double N) {
    Eigen::VectorXd v(size);
    for (size_t i=0; i<size; i++) {
        v[i] = 0;
        double a = basis.supp_min(i);
        double b = basis.supp_max(i);
        for (size_t j=0; j<10; j++) {
            double e = (b-a)/2 *gaussX_10[j] + (a+b)/2;                 
            v[i] +=  gaussW_10[j]*exp(-e/T)*basis(i, e)*pow(e,0.5);; 
        }
        v[i] *= (b-a)/2;
        v[i] *= N*pow(T, -1.5)*2/pow(Constant::Pi,0.5);
    }
    Eigen::VectorXd u = this->basis.Sinv(v);
    for (size_t i=0; i<size; i++) {
        this->f[i] += u[i];
    }
}


// void Distribution::add_density_distribution(vector<double> densities){
//     assert(densities.size() == size);
//     Eigen::VectorXd v(size);
//     std::vector<double> energies = get_trimmed_knots(get_knot_energies());//basis.avg_e;//get_trimmed_knots(get_knot_energies());
//     for (size_t i=0; i<size; i++) {
//         v[i] = 0;
//         double a = basis.supp_min(i);
//         double b = basis.supp_max(i);
//         v[i] += densities[i];//*basis(i, energies[i]);
//         v[i] *= basis.areas[i];//*=(b-a)/2;
//     }
//     Eigen::VectorXd u = this->basis.Sinv(v);
//     for (size_t i=0; i<size; i++) {
//         this->f[i] += u[i];
//     }
// }

// Assumes new basis has same B spline order.
void Distribution::transform_basis(std::vector<double> new_knots){
    int new_basis_order = basis.BSPLINE_ORDER;
    //// Get knots that have densities
    int num_new_splines = get_trimmed_knots(new_knots).size();   // TODO replace get_trimmed_knots with get_num_funcs?. 
    //// Compute densities for knots
    std::vector<std::vector<double>> new_densities(num_new_splines, std::vector<double>(64, 0));

    // Cackle and iterate through each new spline.
    for (size_t i=0; i<num_new_splines; i++){
        // Use current basis to generate the density terms for gaussian integration at for each basis point.   
        // Black magic. ଘ(੭ˊᵕˋ)੭.*･｡ﾟ
        double a = new_knots[i];                  // i.e. <new_basis>.supp_min(i);
        double b = new_knots[i+new_basis_order];  // i.e. <new_basis>.supp_max(i);        
        for(size_t j=0; j < 64; j++){
            double e = (b-a)/2 *gaussX_64[j] + (a+b)/2;
            new_densities[i][j] = (*this)(e);  
        }
    }
    // Change distribution to empty one in new basis.    
    vector<double> new_f(num_new_splines,0);
    set_distribution(new_knots,new_f);
    // Add densities in new basis.
    add_density_distribution(new_densities);
}

void Distribution::add_density_distribution(vector<vector<double>> densities){

    assert(densities.size() == size);
    Eigen::VectorXd v(size);
    std::vector<double> energies = get_trimmed_knots(get_knot_energies());//basis.avg_e;//get_trimmed_knots(get_knot_energies());
    for (size_t i=0; i<size; i++) {
        v[i] = 0;
        double a = basis.supp_min(i);
        double b = basis.supp_max(i);
        // This uses gaussian quadrature (https://pomax.github.io/bezierinfo/legendre-gauss.html)
        // to approximate the integral of bspline_i.density_i between the boundaries of the splines, but with some basis shenanigans to spice it up.
        // Thus v[i] would correspond to the total num electrons contributed by spline i. -S.P.
        for (size_t j=0; j<64; j++) {
            double e = (b-a)/2 *gaussX_64[j] + (a+b)/2;    // e = element E_i of gaussian quadrature sum
            v[i] += gaussW_64[j]*basis.raw_bspline(i, e)*densities[i][j]; 
        }
        v[i] *= (b-a)/2;//*=densities[i]*(b-a)/2;
    }
    // Su = v, where u[i] is the electron density, S is the (sparse) overlap matrix, and v[i] is the num electrons the spline represents. 
    Eigen::VectorXd u = this->basis.Sinv(v);
    for (size_t i=0; i<size; i++) {
        this->f[i] += u[i];
    }
}

std::vector<double> Distribution::get_trimmed_knots(std::vector<double> knots){
    // Remove boundary knots
    while(knots[0] <= basis.min_elec_e() && knots.size() > 0){
        knots.erase(knots.begin());
    }        
    while (knots.back() > basis.max_elec_e() && knots.size() > 0){
        knots.erase(knots.end() - 1);
    }   
    return knots;
}

// ==========================
// IO

// Returns knot energies in eV
std::string Distribution::output_knots_eV() {
    std::stringstream ss;
    for (size_t i=0;i<basis.gridlen(); i++) {
        ss << basis.grid(i)*Constant::eV_per_Ha<<" ";
    }
    return ss.str();
}

// Returns energies in eV
std::string Distribution::output_energies_eV(size_t num_pts) {
    std::stringstream ss;
    size_t pts_per_knot = num_pts / basis.gridlen();
    if (pts_per_knot == 0){
        pts_per_knot = 1;
        cerr<<"[ warn ] Outputting a small number of points per knot."<<endl;
    }
    for (size_t i=0; i<basis.gridlen()-1; i++){
        double e = basis.grid(i);
        double de = (basis.grid(i+1) - e)/(pts_per_knot-1);
        for (size_t j=0; j<pts_per_knot; j++) {
            ss << e*Constant::eV_per_Ha<<" ";
            e += de;
        }
    }
    
    return ss.str();
}

std::string Distribution::output_densities(size_t num_pts) const {
    std::stringstream ss;
    const double units = 1./Constant::eV_per_Ha/Constant::Angs_per_au/Constant::Angs_per_au/Constant::Angs_per_au;
    size_t pts_per_knot = num_pts / basis.gridlen();
    if (pts_per_knot == 0){
        pts_per_knot = 1;
        cerr<<"[ warn ] Outputting a small number of points per knot."<<endl;
    }
    for (size_t i=0; i<basis.gridlen()-1; i++){
        double e = basis.grid(i);
        double de = (basis.grid(i+1) - e)/(pts_per_knot-1);
        for (size_t j=0; j<pts_per_knot; j++) {
            ss << (*this)(e)*units<<" ";
        e += de;
        }
    }
    return ss.str();
}


/**
 * @brief Returns the spline-interpolated electron density (in atomic units, Ha/a0^-3 [Hartree]/[bohr radius]^3) at any energy.
 * 
 * @details This does a dot product of the "Frobenius-treated" basis with each spline's expansion coefficient.
 * @param e energy to get density from.
 * @return double 
 */
double Distribution::operator()(double e) const{
    double tmp=0;
    for (size_t j = 0; j < size; j++) {
        tmp += basis(j, e)*f[j];
    }
    return tmp;
}
/*
void Distribution::addDeltaSpike(double e, double N) {
    int idx = basis.lower_i_from_e(e);
    assert(idx >= 0 && (unsigned) idx+1 < size);
    double A1 = basis.areas[idx];
    double A2 = basis.areas[idx+1];
    double E1 = basis.avg_e[idx];
    double E2 = basis.avg_e[idx+1];
    // Solves the matrix equation
    // A1 * a + A2 * b = N
    // A1 E1 a + A2 E2 b = (A1 + A2) e
    
    double det = (E2 - E1)*A1*A2;
    // Inverse matrix is 
    // 1/det * [ A2 E2   -A2 ]
    //         [-A1 E1    A1 ]
    
    f[idx] +=     ( A2 * E2 * N   -  A2 * (A1 + A2)* e) / det;
    f[idx + 1] += (- A1 * E1 * N  +  A1 * (A1 + A2)* e) / det;
}*/

void Distribution::addDeltaSpike(double e, double N) {
    int idx = basis.i_from_e(e);
    f[idx] += N/basis.areas[idx];
}

/**
 * @brief f <-- f + df/dt|basis  
 * 
 * @param v 
 */
void Distribution::applyDelta(const Eigen::VectorXd& v) {
    Eigen::VectorXd u(size);
    u= (this->basis.Sinv(v)); 
    for (size_t i=0; i<size; i++) {
        f[i] += u[i];
    }
}


// - 3/sqrt(2) * 3 sqrt(e) * f(e) / R_
// Very rough approximation used here
void Distribution::addLoss(const Distribution& d, const LossGeometry &l, double rho) {
    // f += "|   i|   ||   |_"

    double escape_e = 1.333333333*Constant::Pi*l.L0*l.L0*rho;
    for (size_t i=basis.i_from_e(escape_e); i<size; i++) {
        
        f[i] -= d[i] * sqrt(basis.avg_e[i]/2) * l.factor();
    }
}

void Distribution::addDeltaLike(Eigen::VectorXd& v, double e, double height) {
    assert(e>basis.supp_min(0));
    for (size_t i=0;i<size; i++) {
        v[i] += basis(i, e) * height;
    }

}

double Distribution::norm() const{
    assert(f.size() == Distribution::size);  // A very blessed check, lord thank you Sanders -S.P.
    double x=0;
    for (auto&fi : f) {
        x+= fabs(fi);
    }
    return x;
}

double Distribution::density() const {
    double tmp=0;
    for (size_t i = 0; i < size; i++)
    {
        tmp += basis.areas[i]*f[i];
    }
    return tmp;
}

double Distribution::density(size_t cutoff) const {
    double tmp=0;
    for (size_t i = 0; i < cutoff; i++)
    {
        tmp += basis.areas[i]*f[i];
    }
    return tmp;
}

// Returns an estimate of k*T
// based on kinetic energy density of the plasma below the index given by 'cutoff'
double Distribution::k_temperature(size_t cutoff) const {
    double tmp=0;
    // Dodgy integral of e * f(e) de
    double n = this->density(cutoff);
    for (size_t i = 0; i < cutoff; i++)
    {
        tmp += basis.avg_e[i]*f[i]*basis.areas[i];
    }
    return tmp*2./3./n;
}

// Returns an estimate of ln(N_D) based on density and low-energy arguments
double Distribution::CoulombLogarithm() const {
    double tmp=0;
    // Dodgy integral of e * f(e) de
    double n = 0;
    for (size_t i = 0; i < CoulombLog_cutoff; i++) {
        tmp += basis.avg_e[i]*f[i]*basis.areas[i];
        n += f[i]*basis.areas[i];
    }
    double kT = tmp*2./3./n;
    double DebyeLength3 = pow(kT/4/Constant::Pi/n,1.5);
    n = this->density();
    if (n < this->CoulombDens_min) return 0;
    return log(4.*Constant::Pi* n*DebyeLength3);
}


double Distribution::integral(double (g)(double)) {
    double retval = 0;
    for (size_t J = 0; J < Distribution::size; J++)
    {
        double a = basis.supp_min(J);
        double b = basis.supp_max(J);
        double tmp = 0;
        for (size_t i=0;i<10; i++)
        {
            double e = gaussX_10[i]*(b-a)/2 + (b+a)/2;
            tmp += gaussW_10[i]*basis(J, e)*g(e);
        }
        retval +=  tmp * f[J] * (b-a)/2;
    }
    return retval;
}

ostream& operator<<(ostream& os, const Distribution& dist) {
    for (size_t i= 0; i < Distribution::size; i++) {
        os<<" "<<dist[i]/Constant::eV_per_Ha;
    }
    return os;
}

