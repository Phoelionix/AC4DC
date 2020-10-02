#include "FreeDistribution.h"
#include "Dipole.h"
#include "Constant.h"
#include "SplineIntegral.h"

// #define NDEBUG
// to remove asserts

// Initialise static things
size_t Distribution::size=0;
SplineIntegral Distribution::basis;

// Psuedo-constructor thing
void Distribution::set_elec_points(size_t n, double min_e, double max_e, GridSpacing grid_style) {
    // Defines a grid of n points
    basis.set_parameters(n, min_e, max_e, grid_style);
    Distribution::size=n;
}

// Adds Q_eii to the parent Distribution
void Distribution::get_Q_eii (Eigen::VectorXd& v, size_t a, const bound_t& P) const {
    assert(basis.has_Qeii());
    assert(P.size() == basis.Q_EII[a].size());
    assert((unsigned) v.size() == size);
    
    for (size_t xi=0; xi<P.size(); xi++) {
        // Loop over configurations that P refers to
        for (size_t J=0; J<size; J++) {
            for (size_t K=0; K<size; K++) {
                v[J] += P[xi]*f[K]*basis.Q_EII[a][xi][J][K];
            }
        }
    }
}

// Puts the Q_TBR changes in the supplied vector v
void Distribution::get_Q_tbr (Eigen::VectorXd& v, size_t a, const bound_t& P) const {
    assert(basis.has_Qtbr());
    assert(P.size() == basis.Q_TBR[a].size());
    for (size_t eta=0; eta<P.size(); eta++) {
        // Loop over configurations that P refers to
        for (size_t J=0; J<size; J++) {
            for (auto& q : basis.Q_TBR[a][eta][J]) {
                 v[J] += q.val * P[eta] * f[q.K] * f[q.L];
            }
        }
    }
}


// Puts the Q_EE changes into v
void Distribution::get_Q_ee(Eigen::VectorXd& v) const {
    assert(basis.has_Qee());
    // KLUDGE: Fix DebyeLength at 5 Angstrom = 9.4 Bohr
    // const double DebyeLength = 5. / Constant::Angs_per_au;
    // double CoulombLog = log(4./3.*Constant::Pi*DebyeLength*DebyeLength*DebyeLength*density());
    double CoulombLog=5;
    // double CoulombLog = CoulombLogarithm(size/3);
    if (isnan(CoulombLog) || CoulombLog <= 0) return;
    for (size_t J=0; J<size; J++) {
        for (size_t K=0; K<size; K++) {
            for (auto& q : basis.Q_EE[J][K]) {
                 v[J] += q.val * f[K] * f[q.idx] * CoulombLog ;
            }
        }
    }
}

// // Taken verbatim from Rockwood as quoted by Morgan and Penetrante in ELENDIF
// void Distribution::add_Q_ee(const Distribution& d, double kT) {
//     double density=0;
//     double CoulombLog = log(kT/(4*Constant::Pi*density));
//     double alpha = 2*Constant::Pi*sqrt(2)/3*CoulombLog;


// }

// ============================
// utility

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
    double e = basis.min_elec_e();
    double de = (basis.max_elec_e() - e)/(num_pts-1);
    for (size_t i=0;i<num_pts; i++) {
        ss << e*Constant::eV_per_Ha<<" ";
        e += de;
    }
    return ss.str();
}

std::string Distribution::output_densities(size_t num_pts) const {
    std::stringstream ss;
    double e = basis.min_elec_e();
    double de = (basis.max_elec_e() - e)/(num_pts-1);
    for (size_t i=0;i<num_pts; i++) {
        ss << (*this)(e)/Constant::eV_per_Ha<<" ";
        e += de;
    }
    return ss.str();
}

double Distribution::operator()(double e) const{
    double tmp=0;
    for (size_t j = 0; j < size; j++) {
        tmp += basis(j, e)*f[j];
    }
    return tmp;
}

void Distribution::addDeltaSpike(double e, double N) {
    int idx = basis.i_from_e(e);
    // assert(idx >= 0 && (unsigned) idx < size);
    // f[idx] += N/basis.areas[idx];
    unsigned num_offset = basis.BSPLINE_ORDER;
    assert(idx + num_offset < size && idx >= num_offset);
    f[idx] += N/basis.areas[idx];
    // f[idx-1] += N/basis.areas[idx-1]/4;
    // f[idx+1] += N/basis.areas[idx+1]/4;
}

void Distribution::applyDelta(const Eigen::VectorXd& v) {
    Eigen::VectorXd u(size);
    u= (this->basis.Sinv(v));
    for (size_t i=0; i<size; i++) {
        f[i] += u[i];
    }
}


// - 3/sqrt(2) * 3 sqrt(e) * f(e) / R_
// Very rough approximation used here
void Distribution::addLoss(const Distribution& d, const LossGeometry &l) {
    // f += "|   i|   ||   |_"
    for (size_t i=0; i<size; i++) {
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
    assert(f.size() == Distribution::size);
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
double Distribution::CoulombLogarithm(size_t cutoff) const {
    double tmp=0;
    // Dodgy integral of e * f(e) de
    double n = this->density(cutoff);
    for (size_t i = 0; i < cutoff; i++)
    {
        tmp += basis.avg_e[i]*f[i]*basis.areas[i];
    }
    double kT = tmp*2./3./n;
    double DebyeLength3 = pow(kT/4/Constant::Pi/n,1.5);
    return log(4./3.*Constant::Pi* n*DebyeLength3);
    // return 10;
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
