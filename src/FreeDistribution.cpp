#include "FreeDistribution.h"
#include "Dipole.h"
#include "Constant.h"
#include "RateSystem.h"
#include <math.h>
#include <assert.h>

// #define NDEBUG
// to remove asserts

// Initialise static things
size_t Distribution::size=0;
BasisSet Distribution::basis;
Distribution::Q_eii_t Distribution::Q_EII;
Distribution::Q_tbr_t Distribution::Q_TBR;

// Psuedo-constructor thing
void Distribution::set_elec_points(size_t n, double min_e, double max_e){
    // Defines a grid of n points
    basis.set_parameters(n, min_e, max_e);
    Distribution::size=n;
}

// Adds Q_eii to the parent Distribution
void Distribution::apply_Q_eii (Eigen::VectorXd& v, size_t a, const bound_t& P) const {
    for (int xi=0; xi<P.size(); xi++){
        // Loop over configurations that P refers to
        for (int J=0; J<size; J++){
            for (int K=0; K<size; K++){
                v[J] += P[xi]*f[K]*Q_EII[a][xi][J][K];
            }
        }
    }
}

// Adds Q_tbr to the parent Distribution
void Distribution::apply_Q_tbr (Eigen::VectorXd& v, size_t a, const bound_t& P) const {
    for (int xi=0; xi<P.size(); xi++){
        // Loop over configurations that P refers to
        for (int J=0; J<size; J++){
            for (int L=0; L<size; L++){
                for (int K=0; K<size; K++){
                    v[J] += P[xi]*f[L]*f[K]*Q_TBR[a][xi][L][J][K];
                }
            }
        }
    }
}

// Taken verbatim from Rockwood as quoted by Morgan and Penetrante in ELENDIF
// Stores in f the value of Qee[f] (atomic units)
void Distribution::apply_Qee(Eigen::VectorXd& v) const {
    // for (int J=0; J<size; J++){
    //     for (int L=0; L<size; L++){
    //         for (int K=0; K<size; K++){
    //             v[J] += Q_EE[a][xi][J][L][K]*f[K]*f[K];
    //         }
    //     }
    // }
}

//============================
// utility

void Distribution::set_maxwellian(double N, double T){
    T *= Constant::kb_eV; // convert T to units of eV
    Eigen::VectorXd v(size);
    for (size_t i=0; i<size; i++){
        f[i] = 0;
        double a = basis.supp_min(i);
        double b = basis.supp_max(i);
        for (int j=0; j<10; j++){
            double e = (b-a)/2 *gaussX_10[j] + (a+b)/2;
            v[i] += (b-a)/2 *gaussW_10[j]*exp(-e/T)*basis(i, e)/pow(e*Constant::Pi*T, 0.5);
        }
    }
    v = this->basis.Sinv(v);
    for (size_t i=0; i<size; i++){
        f[i] = v[i]*N;
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


// size_t Distribution::i_from_e(double e){
//     // performs a binary search for the basis func with the biggest peak
//     long i = 0;
//     if (i < 0) return 0;
//     if (i >= size) return size-i;
//     return i;
// }

void Distribution::applyDelta(const Eigen::VectorXd& v){
    Eigen::VectorXd u(size);
    u= (this->basis.Sinv(v));
    for (size_t i=0; i<size; i++){
        f[i] += u[i];
    }
}

void Distribution::addDeltaLike(Eigen::VectorXd& v, double e, double height){
    assert(e>basis.supp_min(0));

    for (int i=0; i<size; i++){
        v[i] += basis(i, e) * height;
    }

}

double Distribution::norm() const{
    double x=0;
    for (auto&fi : f){
        x+= fabs(fi);
    }
    return x;
}


// ====================================================================
// Precalculators


// Resizes containers
void Distribution::precompute_Q_coeffs(vector<RateData::Atom>& Store){
    Q_EII.resize(state_type::num_atoms());
    Q_TBR.resize(state_type::num_atoms());
    for (size_t a=0; a<state_type::num_atoms(); a++){
        Q_EII[a].resize(state_type::P_size(a));
        for (auto& Qa : Q_EII[a]){
            Qa.resize(size);
            for (auto& QaJ : Qa){
                    QaJ.resize(size);
            }
        }
        // NOTE: length of EII vector is generally shorter than length of P
        // (usually by 1 or 2, so dense matrices are fine)
        for (auto& eii : Store[a].EIIparams){
            for (size_t J=0; J<size; J++){
                for (size_t K=0; K<size; K++){
                    Q_EII[a][eii.init][J][K] = calc_Q_eii(eii, J, K);
                }
            }
        }
    }
}

void Distribution::Gamma_eii(Eigen::SparseMatrix<double>& Gamma, const RateData::EIIdata& eii, size_t K) {
    // 4pi*sqrt(2) \int_0^\inf e^1/2 phi_k(e) sigma(e) de
    for (size_t eta = 0; eta < eii.fin.size(); eta++) {
        double a = basis.supp_min(K);
        double b = basis.supp_max(K);
        double tmp=0;
        for (int i=0; i<10; i++){
            double e = gaussX_10[i]*(b-a)/2 + (a+b)/2;
            tmp += gaussW_10[i]* basis(K, e)*pow(e,0.5)*Dipole::sigmaBEB(e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
        }
        tmp *= (b-a)/2;
        tmp *= 4*Constant::Pi*1.4142; // Electron mass = 1 in atomic units
        Gamma.coeffRef(eii.fin[eta], eii.init) += tmp;
    }
}

void Distribution::Gamma_tbr(Eigen::SparseMatrix<double>& Gamma, const RateData::EIIdata& eii, size_t J, size_t K) {
    for (size_t eta = 0; eta<eii.fin.size(); eta++) {
        double tmp=0;
        double B=eii.ionB[eta];

        double aK = basis.supp_min(K);
        double bK = basis.supp_max(K);
        double aJ = basis.supp_min(J);
        double bJ = basis.supp_max(J);

        for (int i=0; i<10; i++){
            double e = gaussX_10[i]*(bK-aK)/2 + (aK+bK)/2;
            double tmp2=0;

            for (int j=0; j<10; j++){
                double s = gaussX_10[j]*(bJ-aJ)/2 + (aJ+bJ)/2;
                double ep = s + e + B;
                tmp2 += gaussW_10[j]*basis(J, s)*ep*pow(s,-0.5)*
                    Dipole::DsigmaBEB(ep, e, B, eii.kin[eta], eii.occ[eta]);
            }
            tmp2 *= (bJ-aJ)/2;
            tmp += gaussW_10[i]*basis(K, e)*tmp2/e;
        }
        tmp *= pow(2*Constant::Pi, 4) / 1.4142;
        tmp *= (bK-aK)/2;

        Gamma.coeffRef(eii.init, eii.fin[eta]) += tmp;
    }
}

double Distribution::calc_Q_eii( const RateData::EIIdata& eii, size_t J, size_t K) {
    // computes sum of all Q_eii integrals away from eii.init, integrated against basis functions f_J, f_K
    // J is the 'free' index
    // sqrt(2) sum_eta \int dep sqrt(ep) f(ep) dsigma/de (e | ep) - sqrt(e) sigma * f_J(e)
    // 'missing' factor of 2 in front of RHS is not a mistake! (arises from definition of dsigma)
    double aK = basis.supp_min(K);
    double bK = basis.supp_max(K);
    double aJ = basis.supp_min(J);
    double bJ = basis.supp_max(J);

    double retval=0;
    // Sum over all transitions to eta values
    for (size_t eta = 0; eta<eii.fin.size(); eta++) {
        for (int j=0; j<10; j++){
            double e = gaussX_10[j]*(bJ-aJ)/2 + (aJ+bJ)/2;
            double tmp=0;
            for (int k=0; k<10; k++){
                double ep = gaussX_10[k]*(bK-aK)/2 + (aK+bK)/2;
                tmp += gaussW_10[k]*basis(K, ep)*pow(ep,0.5)*
                    Dipole::DsigmaBEB(ep, e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
            }
            retval += gaussW_10[j]*basis(J, e)*tmp;
        }
        // RHS part
        retval *= (bJ-aJ)/2;
        retval *= (bK-aK)/2;
        double tmp=0;
        for (int k=0; k<10; k++){
            double e = gaussX_10[k]*(bK-aK)/2 + (aK+bK)/2;
            tmp += pow(e, 0.5)*basis(K, e)*Dipole::sigmaBEB(e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
        }
        tmp *= (bK-aK)/2;
        retval -= tmp;
    }

    retval *= 1.4142135624;

    return retval;
}

ostream& operator<<(ostream& os, const Distribution& dist){
    for (size_t i= 0; i < Distribution::size; i++) {
        os<<" "<<dist[i];
    }
    return os;
}
