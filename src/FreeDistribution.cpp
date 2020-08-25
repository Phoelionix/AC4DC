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
bool Distribution::has_Qeii=false;
bool Distribution::has_Qtbr=false;

// Psuedo-constructor thing
void Distribution::set_elec_points(size_t n, double min_e, double max_e){
    // Defines a grid of n points
    basis.set_parameters(n, min_e, max_e);
    Distribution::size=n;
}

// Adds Q_eii to the parent Distribution
void Distribution::get_Q_eii (Eigen::VectorXd& v, size_t a, const bound_t& P) const {
    assert(has_Qeii);
    assert(P.size() == Q_EII[a].size());
    for (int xi=0; xi<P.size(); xi++){
        // Loop over configurations that P refers to
        for (int J=0; J<size; J++){
            for (int K=0; K<size; K++){
                v[J] += P[xi]*f[K]*Q_EII[a][xi][J][K]*1e-6;
            }
        }
    }
}

// Adds Q_tbr to the parent Distribution
void Distribution::get_Q_tbr (Eigen::VectorXd& v, size_t a, const bound_t& P) const {
    assert(has_Qtbr);
    for (int eta=0; eta<P.size(); eta++){
        // Loop over configurations that P refers to
        for (int J=0; J<size; J++){
            for (auto& q : Q_TBR[a][eta][J]){
                 v[J] += q.val * P[eta] * f[q.K] * f[q.L];
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

// Expects T to have units of Kelvin
void Distribution::set_maxwellian(double N, double T){
    T *= Constant::kb_Ha;
    Eigen::VectorXd v(size);
    for (size_t i=0; i<size; i++){
        v[i] = 0;
        double a = basis.supp_min(i);
        double b = basis.supp_max(i);
        for (int j=0; j<10; j++){
            double e = (b-a)/2 *gaussX_10[j] + (a+b)/2;
            double x = N*gaussW_10[j]*exp(-e/T)*basis(i, e)*4*Constant::Pi*pow(2*e, 0.5)*pow(2*Constant::Pi*T,-1.5);
            v[i] += (b-a)/2 * x;
        }
    }
    Eigen::VectorXd u = this->basis.Sinv(v);
    for (size_t i=0; i<size; i++){
        this->f[i] = u[i];
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
    assert(f.size() == Distribution::size);
    double x=0;
    for (auto&fi : f){
        x+= fabs(fi);
    }
    return x;
}


// ====================================================================
// Precalculators


// Resizes containers and fills them with the appropriate values
void Distribution::precompute_Q_coeffs(vector<RateData::Atom>& Store){
    std::cout<<"[ Q precalc ] Beginning coefficient computation..."<<std::endl;
    Q_EII.resize(state_type::num_atoms());
    Q_TBR.resize(state_type::num_atoms());
    for (size_t a=0; a<state_type::num_atoms(); a++){
        std::cout<<"[ Q precalc ] Atom "<<a+1<<"/"<<Store.size()<<std::endl;
        Q_EII[a].resize(state_type::P_size(a));
        Q_TBR[a].resize(state_type::P_size(a));
        for (size_t eta=0; eta<state_type::P_size(a); eta++){
            Q_EII[a][eta].resize(size);
            Q_TBR[a][eta].resize(size);
            for (auto& QaJ : Q_EII[a][eta]){
                    QaJ.resize(size, 0);
            }
            for (auto& QaJ : Q_TBR[a][eta]){
                    QaJ.resize(0); // thse are the sparse lists
            }
        }
        // NOTE: length of EII vector is generally shorter than length of P
        // (usually by 1 or 2, so dense matrices are fine)

        #pragma omp parallel
        #pragma omp for
        for (size_t J=0; J<size; J++){
            std::cout<<"[ Q precalc ] thread "<<omp_get_thread_num()<<std::endl;
            for (size_t K=0; K<size; K++){
                for (auto& eii : Store[a].EIIparams){
                    Q_EII[a][eii.init][J][K] = calc_Q_eii(eii, J, K);
                }
            }
            for (auto& tbr : RateData::inverse(Store[a].EIIparams)){
                Q_TBR[a][tbr.fin][J] = calc_Q_tbr(tbr, J);
            }
        }
    }
    #ifdef DEBUG
    for (size_t a=0; a<state_type::num_atoms(); a++){
        std::cerr<<"[ DEBUG ] Atom "<<a<<std::endl;
        for (size_t i=0; i<state_type::P_size(a); i++){
            std::cerr<<"[ DEBUG ] Config "<<i<<std::endl;
            for (size_t J=0; J<size; J++){
                std::cerr<<"[ DEBUG ] eii>  ";
                for (size_t K=0; K<size; K++){
                    std::cerr<<Q_EII[a][i][J][K]<<" ";
                }
                std::cerr<<std::endl;
                std::cerr<<"[ DEBUG ] tbr>  ";
                for (auto& tup : Q_TBR[a][i][J]){
                    std::cerr<<"("<<tup.K<<","<<tup.L<<","<<tup.val<<") ";
                }
                std::cerr<<std::endl;
            }
        }
    }
    #endif
    std::cout<<"[ Q precalc ] Done."<<std::endl;
    has_Qeii = true;
}

void Distribution::Gamma_eii(GammaType::eiiGraph& Gamma, const std::vector<RateData::EIIdata>& eiiVec, size_t K) {
    // 4pi*sqrt(2) \int_0^\inf e^1/2 phi_k(e) sigma(e) de
    Gamma.resize(0);
    for (size_t init=0; init<eiiVec.size(); init++){
        assert(eiiVec[init].init == init);
        const RateData::EIIdata& eii = eiiVec[init];
        // index on this vector refers to final states
        std::vector<GammaType::NamedValue> gamma_vec(0);
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
            GammaType::NamedValue rate(eii.fin[eta], tmp);
            gamma_vec.push_back(rate);
        }
        Gamma.push_back(gamma_vec);
    }
}

void Distribution::Gamma_tbr(GammaType::eiiGraph& Gamma, const std::vector<RateData::EIIdata>& eiiVec, size_t J, size_t K) {
    Gamma.resize(0);
    for (size_t init=0; init<eiiVec.size(); init++){
        assert(eiiVec[init].init == init);
        const RateData::EIIdata& eii = eiiVec[init];
        // index on this vector refers to final states
        std::vector<GammaType::NamedValue> gamma_vec(0);
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

            GammaType::NamedValue rate(eii.fin[eta], tmp);
            gamma_vec.push_back(rate);
        }
        Gamma.push_back(gamma_vec);
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
        // RHS part, 'missing' factor of 2 is incorporated into definition ofsigma
        retval *= (bJ-aJ)/2;
        retval *= (bK-aK)/2;
        double tmp=0;
        for (int k=0; k<10; k++){
            double e = gaussX_10[k]*(bK-aK)/2 + (aK+bK)/2;
            tmp += pow(e, 0.5)*basis(J, e)*basis(K, e)*Dipole::sigmaBEB(e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
        }
        tmp *= (bK-aK)/2;
        retval -= tmp;
    }

    retval *= 1.4142135624;
    return retval;
}

// JKL indices, but the Q_TBR coeffs are sparse btw J and L
// (overlap by at most BasisSet::BSPLINE_ORDER either way)
// Returns a vector of pairs (Q-idx, nonzzero-idx)
Distribution::tbr_sparse Distribution::calc_Q_tbr( const RateData::InverseEIIdata& tbr, size_t J) {
    // J is the 'free' index, K and L label the given indices of the free distribution
    //
    double aJ = basis.supp_min(J);
    double bJ = basis.supp_max(J);

    size_t num_nonzero=0;

    Distribution::tbr_sparse retval; // a vector of SparseEntry
    SparseEntry tup;
    for (size_t K=0; K<size; K++){
        for (size_t L=0; L<size; L++){
            double aK = basis.supp_min(K);
            double bK = basis.supp_max(K);
            double aL = basis.supp_min(L);
            double bL = basis.supp_max(L);

            double tmp=0;
            // Sum over all transitions in tbr
            for (size_t xi=0; xi<tbr.init.size(); xi++){
                // Heaviside step function
                double B =tbr.ionB[xi];
                double aJ_eff = (aJ < B) ? B : aJ;
                // First half of integral
                for (int j=0; j<10; j++){
                    double e = gaussX_10[j]*(bJ-aJ_eff)/2 + (aJ_eff+bJ)/2;
                    double tmp2=0;

                    double a = max(aK, e - B - bL); // integral lower bound
                    double b = min(bK, e - B - aL); // integral upper bound
                    if (a >= b) continue;
                    for (int k=0; k<10; k++){
                        double ep = gaussX_10[k]*(b-a)/2 + (b + a)/2;
                        tmp2 += gaussW_10[k]*basis(K,ep)*basis(L,e-ep-B)*Dipole::DsigmaBEB(e, ep, B, tbr.kin[xi], tbr.occ[xi])*pow(e-ep-B,-0.5)/ep;
                    }
                    tmp += gaussW_10[j]*tmp2*pow(e,1.5)*0.5*basis(J,e);
                }
                // Second half of integral
                // f(e) integral0->inf ds
                for (int j=0; j<10; j++){
                    double e = gaussX_10[j]*(bJ-aJ)/2 + (aJ+bJ)/2;
                    double tmp2=0;
                    for (int k=0; k<10; k++){
                        double s = gaussX_10[k]*(bK-aK)/2 + (bK + aK)/2;
                        tmp2 += gaussW_10[k]*(s+e+B)*pow(s,-0.5)*basis(K,s)*Dipole::DsigmaBEB(s+e+B, e, B, tbr.kin[xi], tbr.occ[xi]);
                    }
                    tmp -= gaussW_10[j]*tmp2*pow(e,-0.5)*basis(J,e)*basis(L,e);
                }
            }

            tmp *= 2*Constant::Pi*Constant::Pi;


            if (fabs(tmp) > DBL_CUTOFF_TBR){
                tup.K = K;
                tup.L = L;
                tup.val=tmp;
                retval.push_back(tup);
                num_nonzero++;
            }
        }
    }
    has_Qtbr = true;
    return retval;
}

ostream& operator<<(ostream& os, const Distribution& dist){
    for (size_t i= 0; i < Distribution::size; i++) {
        os<<" "<<dist[i];
    }
    return os;
}
