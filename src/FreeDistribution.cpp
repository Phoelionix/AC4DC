#include "FreeDistribution.h"
#include "Dipole.h"
#include "Constant.h"
#include <math.h>
#include <assert.h>

// #define NDEBUG
// to remove asserts

// Initialise static things
size_t Distribution::size=0;
BasisSet Distribution::basis;


namespace BSpline{
    // Template voodoo stolen from stackexchange
    template <unsigned k>
    double BSpline(double x, double *t)
    {
        if (*t <= x && x < *(t+k))
        {
            double a = (x - *t) / (*(t+k-1) - *t);
            double b = (*(t+k) - x) / (*(t+k) - *(t+1));

            return a * BSpline<k-1>(x, t) + b * BSpline<k-1>(x, (t+1));
        }
        else
            return 0;
    }

    template <>
    double BSpline<1>(double x, double *t)
    {
        if (*t <= x && x < *(t+1))
            return 1.;
        else
            return 0.;
    }

};

// Psuedo-constructor thing
void Distribution::set_elec_points(size_t n, double min_e, double max_e){
    // Defines a grid of n+1 points
    basis.set_parameters(n, min_e, max_e);
    Distribution::size=n;
}

// Computes all EII rates for specified transition for basis func k
// Returns the relevant Gamma rate matrix for atom a
// dP/dt = f*(WP - WT P)
Eigen::SparseMatrix<double> Distribution::Gamma_eii(const RateData::EIIdata& eii, size_t K, size_t a) {
    // 4pi*sqrt(2) \int_0^\inf e^1/2 phi_k(e) sigma(e) de
    Eigen::SparseMatrix<double> Gamma;
    for (size_t xi = 0; xi < eii.fin.size(); xi++) {
        double a = basis.supp_min(K);
        double b = basis.supp_max(K);
        double tmp=0;
        for (int i=0; i<GAUSS_ORDER; i++){
            double e = gaussX_10[i]*(b-a)/2 + (a+b)/2;
            tmp += gaussW_10[i]* basis(K, e)*pow(e,0.5)*Dipole::sigmaBEB(e, eii.ionB[xi], eii.kin[xi], eii.occ[xi]);
        }
        tmp *= (b-a)/2;
        tmp *= 4*Constant::Pi*1.4142; // Electron mass = 1 in atomic units
        Gamma.coeffRef(eii.fin[xi], eii.init) += tmp;
    }
    return Gamma;
}

Eigen::SparseMatrix<double> Distribution::Gamma_tbr(const RateData::EIIdata& eii, size_t J, size_t K, size_t a) {
    Eigen::SparseMatrix<double> Gamma;
    for (size_t xi = 0; xi<eii.fin.size(); xi++) {
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
                tmp2 += gaussW_10[j]*basis(J, s)*ep*pow(s,-0.5)*
                    Dipole::DsigmaBEB(ep, e, B, eii.kin[xi], eii.occ[xi]);
            }
            tmp2 *= (b2-a2)/2;
            tmp += gaussW_10[i]*basis(K, e)*tmp2/e;
        }
        tmp *= pow(2*Constant::Pi, 4) / 1.4142;
        tmp *= (b1-a1)/2;

        Gamma.coeffRef(eii.init, eii.fin[xi]) += tmp;
    }
    return Gamma;
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
    Eigen::VectorXd v(size);
    for (size_t i=0; i<size; i++){
        f[i] = 0;
        double a = basis.supp_min(i);
        double b = basis.supp_max(i);
        for (int j=0; j<GAUSS_ORDER; j++){
            double e = (b-a)/2 *gaussX_10[j] + (a+b)/2;
            v[i] += (b-a)/2 *gaussW_10[j]*exp(-e/T)*basis(i, e)/pow(e*Constant::Pi*T, 0.5);
        }
    }
    v = N*(this->basis.Sinv(v));
    for (size_t i=0; i<size; i++){
        f[i] = v[i];
    }
}

// Returns energies in eV
std::string Distribution::get_energies_eV(){
    std::stringstream ss;
    for (int i=0; i<basis.gridlen(); i++){
        ss << basis.grid(i)*Constant::eV_per_Ha<<" ";
    }
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

void Distribution::addDeltaLike(double e, double height){
    assert(e>basis.supp_min(0));
    Eigen::VectorXd v(size);
    for (int i=0; i<size; i++){
        v[i] = basis(i, e);
    }
    v= height*(this->basis.Sinv(v));
    for (size_t i=0; i<size; i++){
        f[i] += v[i];
    }
}

double Distribution::norm() const{
    double x=0;
    for (auto&fi : f){
        x+= fabs(fi);
    }
    return x;
}


void BasisSet::set_parameters(size_t n, double min, double max) {
    num_funcs = n;
    knot.resize(n+BSPLINE_ORDER);
    knot[0] = min;
    double de = 1.*(max-min)/(n+BSPLINE_ORDER+1);
    for(int i=1; i<n+BSPLINE_ORDER; i++){
        knot[i] = knot[i-1] + de;
    }
    // Compute overlap matrix
    Eigen::SparseMatrix<double> S(num_funcs, num_funcs);
    // Eigen::MatrixXd S(num_funcs, num_funcs);
    for (int i=0; i<num_funcs; i++){
        for (int j=i+1; j<num_funcs; j++){
            double tmp=overlap(i, j);
            if (tmp != 0){
                S.insert(i,j) = tmp;
                S.insert(j,i) = tmp;
            }
        }
        S.insert(i,i) = overlap(i,i);
    }
    // Compute the Cholesky decomposition
    cholmachine.compute(S);
    if(cholmachine.info()!=Eigen::Success) {
      std::cerr<<"Cholesky Factorisation of overlap matrix failed!"<<endl;
    }
    // Dense only
    // if(cholmachine.rcond() < 1e-6) {
    //     std::cerr<<"Condition number is "<<1./cholmachine.rcond();
    //     std::cerr<<"S-matrix inversion may be numerically unstable."<<endl;
    // }

}

Eigen::VectorXd BasisSet::Sinv(Eigen::VectorXd deltaf){
    // Solves the linear system S fdot = deltaf
    auto x = cholmachine.solve(deltaf);
    #ifndef NDEBUG
    if(cholmachine.info()!=Eigen::Success) {
      std::cerr<<"Linear system could not be solved!"<<std::endl;
    }
    #endif
    return x;
}

double BasisSet::operator()(size_t i, double x){
    // Returns the i^th B-spline of order BSPLINE_ORDER
    return BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
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
        double e = gaussX_10[i]*(b-a)/2 + (b+a)/2;
        tmp += gaussW_10[i]*(*this)(j, e)*(*this)(k, e);
    }
    tmp *= (b-a)/2;
    return tmp;
}

ostream& operator<<(ostream& os, const Distribution& dist){
    for (size_t i= 0; i < Distribution::size; i++) {
        os<<" "<<dist[i];
    }
    return os;
}
