#include "FreeDistribution.h"
#include "Dipole.h"
#include "Constant.h"
#include <math.h>
// #include <stringstream>
// #define NDEBUG


size_t Distribution::size=0;

// vector<double> Distribution::grid   = vector<double>(0);
// vector<double> Distribution::grid12  = vector<double>(0);
// vector<double> Distribution::widths = vector<double>(0);
// double Distribution::max_e = 2e2f;
// double Distribution::min_e = 10.f;

// This is bad style.
BasisSet* basis = NULL;

// Psuedo-constructor thing
void Distribution::set_elec_points(size_t n, double min_e, double max_e){
    // Defines a grid of n+1 points
    basis = new BasisSet(n, min_e, max_e);
    size=n;
}

// Computes all EII rates for specified transition for basis func k
// Returns the relevant Gamma rate matrix for atom a
// dP/dt = f*(WP - WT P)
void Distribution::Gamma_eii(Eigen::SparseMatrix<double>& Gamma; const CustomDataType::EIIdata& eii, size_t K, size_t a) {
    // 4pi*sqrt(2) \int_0^\inf e^1/2 phi_k(e) sigma(e) de

    for (size_t xi = 0; xi < eii.fin.size(); xi++) {
        double a = basis.supp_min(K);
        double b = basis.supp_max(K);
        double tmp=0;
        for (int i=0; i<GAUSS_ORDER; i++){
            double e = gaussX_10[i]*(b-a)/2 + (a+b)/2
            tmp += gaussW_10[i]* (*basis(K, E))*pow(e,0.5)*Dipole::sigmaBEB(e, eii.ionB[xi], eii.kin[xi], eii.occ[xi]);
        }
        tmp *= (b-a)/2;
        tmp = basis.integral(integrand, K);
        tmp *= 4*Constant::Pi*1.4142; // Electron mass = 1 in atomic units
        Gamma(eii.fin[xi], eii.init) = tmp;
    }
    return Gamma;
}

void Distribution::Gamma_tbr(Eigen::SparseMatrix<double>& Gamma; const CustomDataType::EIIdata& eii, size_t J, size_t K, size_t a) {
    for (size_t xi; xi<eii.fin.size(); xi++) {
        double tmp=0;
        double B=eii.ionB[xi];

        double a1 = basis.supp_min(K);
        double b1 = basis.supp_max(K);
        double a2 = basis.supp_min(J);
        double b2 = basis.supp_max(J);

        for (int i=0; i<GAUSS_ORDER; i++){
            double e = gaussX_10[i]*(b1-a1)/2 + (a1+b1)/2;
            double tmp2=0;

            for (int j=0; j<GAUSS_ORDER; j++){
                double s = gaussX_10[i]*(b2-a2)/2 + (a2+b2)/2;
                double ep = s + e + B;
                tmp2 += gaussW_10[j]*(*basis(J, s))*ep*pow(s,-0.5)*
                    Dipole::DsigmaBEB(ep, e, B, eii.kin[xi], eii.occ[xi]);
            }
            tmp2 *= (b2-a2)/2;
            tmp += gaussW_10[i]*(*basis(K, E))*tmp2/e;
        }
        tmp *= pow(2*Constant::Pi, 4) / 1.4142;
        tmp *= (b1-a1)/2;

        Gamma(ep.init, ep.fin[xi]) = tmp;
    }
}

// Adds Q_eii to the parent Distribution
void Distribution::add_Qeii (size_t a, const Distribution& F, const bound_t& P) {
    // for (int xi=0; xi<P.size(); xi++){
    //     // Loop over configurations that P refers to
    //     for (int i=0; i<size; i++){
    //         this->f[i] += Q_EII[a][xi][i]*P[xi]*F.f[i];
    //     }
    // }
}

// Adds Q_tbr to the parent Distribution
void Distribution::add_Qtbr (size_t a, const Distribution& F, const bound_t& P) {
    // for (int i=0; i<size; i++){
    //     // For each energy value...
    //     for (int xi=0; xi<P.size(); xi++){
    //         // Loop over configurations that P refers to
    //         // for (int j=0; j<size; j++){
    //         //     this->f[i] += P[xi]*Q_TBR[a][xi][i][j]*F.f[i-j-B_idx]*F.f[j]
    //         // }
    //         // n^2 complexity :/
    //     }
    // }
}

// Taken verbatim from Rockwood as quoted by Morgan and Penetrante in ELENDIF
// Stores in f the value of Qee[f] (atomic units)
void Distribution::add_Qee(const Distribution& F) {
    // TBD
}

void Distribution::set_maxwellian(double N, double T){
    T *= Constant::kb_eV; // convert T to units of eV
    for (size_t i=0; i<size; i++){
        f[i] = N*exp(-grid[i]/T)/pow(grid[i]*Constant::Pi*T, 0.5);
    }
}

// Returns energies in eV
std::string Distribution::get_energies_eV(){
    std::stringstream ss;
    // TODO
    return ss.str();
}


//
// double Distribution::e_from_i(size_t i){
//     // Linear grid
//     return min_e + i*(max_e - min_e)/size;
// }
//
// size_t Distribution::i_from_e(double e){
//     long i= std::floor((e - min_e) * size /(max_e - min_e));
//     if (i < 0) return 0;
//     if (i >= size) return size-i;
//     return i;
// }

void Distribution::addDeltaLike(double e, double mass){
    for ( size_t i=0; i<size; i++){
        f[i] += mass*Sinv*basis.val(e);
    }
}


BasisSet::BasisSet(size_t n, double min, double max) : Sinv(n,n), num_funcs(n){
    knot.resize(n+BSPLINE_ORDER);
    knot[0] = min;
    for(int i=1; i<n+BSPLINE_ORDER; i++){
        knot[i] = knot[i-1] + i*(max-min)/(n+BSPLINE_ORDER);
    }
    // Compute overlap matrix
    Eigen::SparseMatrix<double> S;
    for (int i=0; i<size; i++){
        for (int j=i; j<i+BSPLINE_ORDER; j++){
            S.insert(i,j) = overlap(i, j);
            S.insert(j,i) = S.coeffRef(i,j);
        }
    }
    // Compute the Cholesky decomposition
    cholmachine.compute(S);
    if(cholmachine.info()!=Success) {
      std::cerr<<"Cholesky Factorisation of overlap matrix failed!";
    }

}

Eigen::VectorXd BasisSet::Sinv(Eigen::VectorXd deltaf){
    // Solves the linear system S fdot = deltaf
    x = cholmachine.solve(deltaf);
    if(cholmachine.info()!=Success) {
      std::cerr<<"Linear system could not be solved!"
    }
    return x;
}


double BasisSet::operator()(size_t i, double x){
    // Returns the n^th B-spline of order BSPLINE_ORDER
    return BSpline<BSPLINE_ORDER>(x, &knot[i]);
}

inline double BasisSet::supp_max(unsigned k){
    return knot[k+BSPLINE_ORDER];
}

inline double BasisSet::supp_min(unsigned k){
    return knot[k];
}


double BasisSet::overlap(size_t j,size_t k){
    if (j<k){
        int s = k;
        k=j;
        j=s;
    }// j >= k
    if (j-k>BSPLINE_ORDER) return 0;
    double b = supp_max(j);
    double a = supp_min(k);
    double tmp=0;
    for (int i=0; i<GAUSS_ORDER; i++){
        double e = gaussX_10[i](b-a)/2 + (b+a)/2;
        tmp += gaussW_10[i]*(*this)(j, e)*(*this)(k, e);
    }
    tmp *= (b-a)/2;
    return tmp;
}



ostream& operator<<(ostream& os, const Distribution& dist){
    for (size_t i= 0; i < Distribution::size; i++) {
        os<<dist.f[i]<<" ";
    }
    return os;
}
