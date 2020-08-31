#include "SplineIntegral.h"

//                  value to find, sorted array, length of array, index that *arr refers to
int r_binsearch(double val, double* arr, int len, int idx) {
    if (len==1) return idx;
    int mid = len/2;
    // arr has indices i i+1 ... i+len-1
    if (*(arr + mid) > val){
        return r_binsearch(val, arr, len/2, idx);
    } else {
        return r_binsearch(val, arr+len/2, len - (len/2), idx+len/2);
    }
}



int SplineIntegral::i_from_e(double e){
    return r_binsearch(e, this->knot.data(), this->knot.size(), 0);
}

double SplineIntegral::area(size_t J){
    assert(J>=0 && J < num_funcs);
    double min = this->supp_min(J);
    double max = this->supp_max(J);
    double tmp=0;
    for (int i=0; i<10; i++){
        double e = gaussX_10[i]*(max - min)/2 + (max+min)/2;
        tmp += gaussW_10[i]*(*this)(J, e);
    }
    // Numerical integral of a polynomial machine go brrrr
    return tmp*(max-min)/2;
}

// Resizes containers and fills them with the appropriate values
void SplineIntegral::precompute_Q_coeffs(vector<RateData::Atom>& Atoms){
    std::cout<<"[ Q precalc ] Beginning coefficient computation..."<<std::endl;
    Q_EII.resize(Atoms.size());
    Q_TBR.resize(Atoms.size());
    for (size_t a=0; a<Atoms.size(); a++){
        std::cout<<"[ Q precalc ] Atom "<<a+1<<"/"<<Atoms.size()<<std::endl;
        Q_EII[a].resize(Atoms[a].num_conf);
        Q_TBR[a].resize(Atoms[a].num_conf);
        for (size_t eta=0; eta<Atoms[a].num_conf; eta++){
            Q_EII[a][eta].resize(num_funcs);
            Q_TBR[a][eta].resize(num_funcs);
            for (auto& QaJ : Q_EII[a][eta]){
                    QaJ.resize(num_funcs, 0);
            }
            for (auto& QaJ : Q_TBR[a][eta]){
                    QaJ.resize(0); // thse are the sparse lists
            }
        }
        // NOTE: length of EII vector is generally shorter than length of P
        // (usually by 1 or 2, so dense matrices are fine)

        size_t counter=1;
        #pragma omp parallel default(none) shared(a, counter, Atoms, std::cout)
		{
			#pragma omp for schedule(dynamic) nowait
            for (size_t J=0; J<num_funcs; J++){
                #pragma omp critical
                {
                    std::cout<<"[ Q precalc ] "<<counter<<"/"<<num_funcs<<" thread "<<omp_get_thread_num()<<std::endl;
                    counter++;
                }
                for (size_t K=0; K<num_funcs; K++){
                    for (auto& eii : Atoms[a].EIIparams){
                        // XXX: FUDGE FACTOR GOES HERE
                        Q_EII[a][eii.init][J][K] = calc_Q_eii(eii, J, K);
                    }
                }
                for (auto& tbr : RateData::inverse(Atoms[a].EIIparams)){
                    Q_TBR[a][tbr.fin][J] = calc_Q_tbr(tbr, J);
                }
            }
        }
    }
    #ifdef DEBUG
    for (size_t a=0; a<Atoms.size(); a++){
        std::cerr<<"[ DEBUG ] Atom "<<a<<std::endl;
        for (size_t i=0; i<Atoms[a].num_conf; i++){
            std::cerr<<"[ DEBUG ] eii> Config "<<i<<std::endl;
            std::cerr<<"[ DEBUG ] tbr> Config "<<i<<std::endl;
            for (size_t J=0; J<num_funcs; J++){
                std::cerr<<"[ DEBUG ] eii>  ";
                for (size_t K=0; K<num_funcs; K++){
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
    _has_Qeii = true;
    _has_Qtbr = true;
    std::cout<<"[ Q precalc ] Done."<<std::endl;
}

void SplineIntegral::Gamma_eii(eiiGraph& Gamma, const std::vector<RateData::EIIdata>& eiiVec, size_t K) const {
    // 4pi*sqrt(2) \int_0^\inf e^1/2 phi_k(e) sigma(e) de
    Gamma.resize(0);
    for (size_t init=0; init<eiiVec.size(); init++){
        assert(eiiVec[init].init == init);
        const RateData::EIIdata& eii = eiiVec[init];
        // index on this vector refers to final states
        std::vector<SparsePair> gamma_vec(0);
        for (size_t eta = 0; eta < eii.fin.size(); eta++) {
            double a = this->supp_min(K);
            double b = this->supp_max(K);
            double tmp=0;
            for (int i=0; i<10; i++){
                double e = gaussX_10[i]*(b-a)/2 + (a+b)/2;
                tmp += gaussW_10[i]* (*this)(K, e)*pow(e,0.5)*Dipole::sigmaBEB(e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
            }
            tmp *= (b-a)/2;
            tmp *= 4*Constant::Pi*1.4142; // Electron mass = 1 in atomic units
            SparsePair rate(eii.fin[eta], tmp);
            gamma_vec.push_back(rate);
        }
        Gamma.push_back(gamma_vec);
    }
}

void SplineIntegral::Gamma_tbr(eiiGraph& Gamma, const std::vector<RateData::EIIdata>& eiiVec, size_t J, size_t K) const {
    Gamma.resize(0);
    for (size_t init=0; init<eiiVec.size(); init++){
        assert(eiiVec[init].init == init);
        const RateData::EIIdata& eii = eiiVec[init];
        // index on this vector refers to final states
        std::vector<SparsePair> gamma_vec(0);
        for (size_t eta = 0; eta<eii.fin.size(); eta++) {
            double tmp=0;
            double B=eii.ionB[eta];

            double aK = this->supp_min(K);
            double bK = this->supp_max(K);
            double aJ = this->supp_min(J);
            double bJ = this->supp_max(J);

            for (int i=0; i<10; i++){
                double e = gaussX_10[i]*(bK-aK)/2 + (aK+bK)/2;
                double tmp2=0;

                for (int j=0; j<10; j++){
                    double s = gaussX_10[j]*(bJ-aJ)/2 + (aJ+bJ)/2;
                    double ep = s + e + B;
                    tmp2 += gaussW_10[j]*(*this)(J, s)*ep*pow(s,-0.5)*
                        Dipole::DsigmaBEB(ep, e, B, eii.kin[eta], eii.occ[eta]);
                }
                tmp2 *= (bJ-aJ)/2;
                tmp += gaussW_10[i]*(*this)(K, e)*tmp2/e;
            }
            tmp *= pow(2*Constant::Pi, 4) / 1.4142;
            tmp *= (bK-aK)/2;

            SparsePair rate(eii.fin[eta], tmp);
            gamma_vec.push_back(rate);
        }
        Gamma.push_back(gamma_vec);
    }
}

double SplineIntegral::calc_Q_eii( const RateData::EIIdata& eii, size_t J, size_t K) const {
    // computes sum of all Q_eii integrals away from eii.init, integrated against basis functions f_J, f_K
    // J is the 'free' index
    // sqrt(2) sum_eta \int dep sqrt(ep) f(ep) dsigma/de (e | ep) - sqrt(e) sigma * f_J(e)
    // 'missing' factor of 2 in front of RHS is not a mistake! (arises from definition of dsigma)

    double min_J = this->supp_min(J);
    double max_J = this->supp_max(J);

    double retval=0;
    // Sum over all transitions to eta values
    for (size_t eta = 0; eta<eii.fin.size(); eta++) {
        for (int j=0; j<10; j++){
            double e = gaussX_10[j]*(max_J-min_J)/2 + (max_J+min_J)/2;
            double tmp=0;
            double min_ep = max(e+eii.ionB[eta], this->supp_min(K));
            double max_ep = this->supp_max(K);
            if (max_ep <= min_ep) continue;
            for (int k=0; k<10; k++){
                double ep = gaussX_10[k]*(max_ep-min_ep)*0.5 + (min_ep+max_ep)*0.5;
                tmp += gaussW_10[k]*(*this)(K, ep)*pow(ep,0.5)*
                    Dipole::DsigmaBEB(ep, e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
            }
            retval += gaussW_10[j]*(*this)(J, e)*tmp*(max_ep-min_ep);
        }
        // RHS part, 'missing' factor of 2 is incorporated into definition of sigma
        retval *= (max_J-min_J)/2;
        double tmp=0;
        double min_JK = max(this->supp_min(J), this->supp_min(K));
        double max_JK = min(this->supp_max(J), this->supp_max(K));
        if (max_JK <= min_JK) continue;
        for (int k=0; k<10; k++){
            double e = gaussX_10[k]*(max_JK-min_JK)/2 + (min_JK+max_JK)/2;
            tmp += pow(e, 0.5)*(*this)(J, e)*(*this)(K, e)*Dipole::sigmaBEB(e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
        }
        tmp *= (max_JK-min_JK)/2;
        retval -= tmp;
    }

    retval *= 1.4142135624;
    return retval;
}

// JKL indices, but the Q_TBR coeffs are sparse btw J and L
// (overlap by at most BasisSet::BSPLINE_ORDER either way)
// Returns a vector of pairs (Q-idx, nonzzero-idx)
SplineIntegral::sparse_matrix SplineIntegral::calc_Q_tbr( const RateData::InverseEIIdata& tbr, size_t J) const {
    // J is the 'free' index, K and L label the given indices of the free distribution
    //
    double aJ = this->supp_min(J);
    double bJ = this->supp_max(J);

    size_t num_nonzero=0;

    SplineIntegral::sparse_matrix retval; // a vector of SparseTriple
    SparseTriple tup;
    for (size_t K=0; K<num_funcs; K++){
        for (size_t L=0; L<num_funcs; L++){
            double aK = this->supp_min(K);
            double bK = this->supp_max(K);
            double aL = this->supp_min(L);
            double bL = this->supp_max(L);

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
                        tmp2 += gaussW_10[k]*(*this)(K,ep)*(*this)(L,e-ep-B)*Dipole::DsigmaBEB(e, ep, B, tbr.kin[xi], tbr.occ[xi])*pow(e-ep-B,-0.5)/ep;
                    }
                    tmp += gaussW_10[j]*tmp2*pow(e,1.5)*0.5*(*this)(J,e);
                }
                // Second half of integral
                // f(e) integral0->inf ds
                for (int j=0; j<10; j++){
                    double e = gaussX_10[j]*(bJ-aJ)/2 + (aJ+bJ)/2;
                    double tmp2=0;
                    for (int k=0; k<10; k++){
                        double s = gaussX_10[k]*(bK-aK)/2 + (bK + aK)/2;
                        tmp2 += gaussW_10[k]*(s+e+B)*pow(s,-0.5)*(*this)(K,s)*Dipole::DsigmaBEB(s+e+B, e, B, tbr.kin[xi], tbr.occ[xi]);
                    }
                    tmp -= gaussW_10[j]*tmp2*pow(e,-0.5)*(*this)(J,e)*(*this)(L,e);
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
    return retval;
}

SplineIntegral::sparse_matrix SplineIntegral::calc_Q_ee(size_t J) const{
    double min_J = this->supp_min(J);
    double max_J = this->supp_max(J);

    size_t num_nonzero=0;
    double tmp=0;
    SplineIntegral::sparse_matrix retval; // a vector of SparseTriple
    SparseTriple tup;
    for (size_t K=max((size_t) 0, J-BSPLINE_ORDER); K<min(num_funcs, J+BSPLINE_ORDER); K++){
        double min = std::max(this->supp_min(K), min_J);
        double max = std::min(this->supp_max(K), max_J);
        // make sure K and J overlap
        if (fabs(tmp) > DBL_CUTOFF_QEE){
            tup.K = K;
            tup.L = 0;
            tup.val=tmp;
            retval.push_back(tup);
            num_nonzero++;
        }
    }


}
