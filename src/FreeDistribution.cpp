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
void Distribution::set_elec_points(size_t n, double min_e, double max_e, GridSpacing grid_style){
    // Defines a grid of n points
    basis.set_parameters(n, min_e, max_e, 0, grid_style);
    Distribution::size=n;
}

// Adds Q_eii to the parent Distribution
void Distribution::get_Q_eii (Eigen::VectorXd& v, size_t a, const bound_t& P) const {
    assert(basis.has_Qeii());
    assert(P.size() == basis.Q_EII[a].size());
    assert(v.size() == size);
    double dq =0;
    for (int xi=0; xi<P.size(); xi++){
        // Loop over configurations that P refers to
        for (int J=0; J<size; J++){
            for (int K=0; K<size; K++){
                v[J] += P[xi]*f[K]*basis.Q_EII[a][xi][J][K];
            }
        }
    }
}

// Adds Q_tbr to the parent Distribution
void Distribution::get_Q_tbr (Eigen::VectorXd& v, size_t a, const bound_t& P) const {
    assert(basis.has_Qtbr());
    assert(P.size() == basis.Q_TBR[a].size());
    for (int eta=0; eta<P.size(); eta++){
        // Loop over configurations that P refers to
        for (int J=0; J<size; J++){
            for (auto& q : basis.Q_TBR[a][eta][J]){
                 v[J] += q.val * P[eta] * f[q.K] * f[q.L];
            }
        }
    }
}

// Taken verbatim from Rockwood as quoted by Morgan and Penetrante in ELENDIF
// Stores in f the value of Qee[f] (atomic units)
void Distribution::get_Q_ee(Eigen::VectorXd& v) const {
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

// Expects T to have units of Ha
void Distribution::set_maxwellian(double T, double N){
    Eigen::VectorXd v(size);
    for (size_t i=0; i<size; i++){
        v[i] = 0;
        double a = basis.supp_min(i);
        double b = basis.supp_max(i);
        for (int j=0; j<10; j++){
            double e = (b-a)/2 *gaussX_10[j] + (a+b)/2;
            v[i] +=  gaussW_10[j]*exp(-e/T)*basis(i, e)*pow(e,0.5);;
        }
        v[i] *= (b-a)/2;
        v[i] *= N*pow(T, -1.5)*2/pow(Constant::Pi,0.5);
    }
    Eigen::VectorXd u = this->basis.Sinv(v);
    for (size_t i=0; i<size; i++){
        this->f[i] = u[i];
    }
}

// Returns energies in eV
std::string Distribution::output_energies_eV(size_t num_pts) {
    std::stringstream ss;
    double e = basis.min_elec_e();
    double de = (basis.max_elec_e() - e)/(num_pts-1);
    for (int i=0; i<num_pts; i++){
        ss << e*Constant::eV_per_Ha<<" ";
        e += de;
    }
    // for (int i=0; i<basis.gridlen(); i++){
    //     ss << basis.grid(i)*Constant::eV_per_Ha<<" ";
    // }
    return ss.str();
}

std::string Distribution::output_densities(size_t num_pts) const {
    std::stringstream ss;
    double e = basis.min_elec_e();
    double de = (basis.max_elec_e() - e)/(num_pts-1);
    for (int i=0; i<num_pts; i++){
        ss << (*this)(e)/Constant::eV_per_Ha<<" ";
        e += de;
    }
    return ss.str();
}

double Distribution::operator()(double e) const{
    double tmp;
    for (size_t j = 0; j < size; j++) {
        tmp += basis(j, e)*f[j];
    }
    return tmp;
}

void Distribution::addDeltaSpike(double e, double N){
    f[basis.i_from_e(e)] += N/basis.area(basis.i_from_e(e));
}

void Distribution::applyDelta(const Eigen::VectorXd& v){
    Eigen::VectorXd u(size);
    u= (this->basis.Sinv(v));
    for (size_t i=0; i<size; i++){
        f[i] += u[i];
    }
}


// - 3/sqrt(2) * 3 sqrt(e) * f(e) / R_
// Very rough approcimation used here
void Distribution::addLoss(const Distribution d, double R){
    // f += "|   i|   ||   |_"
    for (size_t i=0; i<size; i++){
        f[i] -= d[i] * 3 * sqrt(basis.supp_min(i) + basis.supp_max(i)) /2 / R;
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

double Distribution::integral(double (g)(double)){
    double retval = 0;
    for (size_t J = 0; J < Distribution::size; J++)
    {
        double a = basis.supp_min(J);
        double b = basis.supp_max(J);
        double tmp = 0;
        for (int i=0; i<10; i++)
        {
            double e = gaussX_10[i]*(b-a)/2 + (b+a)/2;
            tmp += gaussW_10[i]*basis(J, e)*g(e);
        }
        retval +=  tmp * f[J] * (b-a)/2;
    }
    return retval;
}


ostream& operator<<(ostream& os, const Distribution& dist){
    for (size_t i= 0; i < Distribution::size; i++) {
        os<<" "<<dist[i]/Constant::eV_per_Ha;
    }
    return os;
}
