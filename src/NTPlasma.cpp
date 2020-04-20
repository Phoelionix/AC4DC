#include "NTPlasma.h"

NTPlasma::NTPlasma(const int T_POINTS, const Grid &lattice): Lattice(lattice) {
    this->e_dist.resize(T_POINTS);
    int N = this->Lattice.size();
    for (size_t i = 0; i < e_dist.size(); i++) {
        this->e_dist[i].resize(N);
    }
}

// df(e,t)/dt = Q_eii[f, p](e) + Q_pi[f](e, t) + Q_a[f](e) + Q_ee[f](e) + Q_tbr[f](e) - Q_loss(e)

// returns integral of g(e) against f
double NTPlasma::fintegral(double (&g)(double))
{
    for (size_t i = 0; i < count; i++) {
        /* code */
    }
}

double eii_aux(double e)
{
    //m_e = 1 in au
    return pow(2*e,1.5)
}

// Electron-impact ionisation collision term at timestep n
double NTPlasma::Q_eii(nte_state_t &f, vector<double> &P, double e){
    // Following Leonov,
    // I_eii = n_a *[ 1/v^2 \sum_l P_l \int v'^3 \sum_mdsigma_lm (v | v') f(v') dv'
    // -v/2 f(v) \sum_l P_l \sum_m sigma_lm ]
    // IN ENERGY COORDS
    // I_eii = n_a *[ m_e/e \sum_l P_l \int (e')^3/2` \sum_m dsigma_lm (e | e') f(e') de'
    // -v/2 f(v) \sum_l P_l \sum_m sigma_lm ]
    double res = 0;
    double tmp;
    for (size_t l = 0; l < P.size(); l++) {
        res += fintegral(eii_aux)*P[l];


        // Integrate over f in the stupid way, since it contains large
        // bumps that even a 13th order GQ integrator may struggle with
    }

}
// Photoionisation flux at timestep n
double NTPlasma::Q_pi(nte_state_t &f, double e);
// Auger effect
double NTPlasma::Q_a(nte_state_t &f, double e);
// Electron-electron scattering
// double Q_ee(int n);
// Density loss due to electrons leaving the system
double NTPlasma::Q_loss(nte_state_t &f, double e);
