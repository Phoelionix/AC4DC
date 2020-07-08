#include "RateSystem.h"
#include "Dipole.h"
#include <math.h>
// #include <stringstream>

vector<double> Distribution::grid   = vector<double>(0);
vector<double> Distribution::widths = vector<double>(0);
vector<size_t> state_type::P_sizes  = vector<size_t>(0);
size_t Distribution::size=0;
double Distribution::max_e = 1e3f;
double Distribution::min_e = 10.f;


void Distribution::set_elec_points(size_t n, double _min_e, double _max_e){
    min_e = _min_e;
    max_e = _max_e;
    size=n;
    grid.resize(n);
    widths.resize(n);

    grid[0] = min_e;
    for (size_t i=1; i<n; i++){
        grid[i] = e_from_i(i);
        widths[i] = grid[i]-grid[i-1];
    }
}

double Distribution::eii_int (const CustomDataType::EIIdata& eii, const int J) const {
    // Returns electrons ionised per atom for dat acorresponding to idx i.
    // TODO: check that this integral makes SENSE
    // INTEGRAL_0^\infty f(e) e/sqrt(2me) SigmaBEB
    double tmp = 0;
    for (size_t i = 0; i < size; i++) {
        tmp += sqrt(grid[i])*f[i]* Dipole::sigmaBEB(grid[i], eii.ionB[J], eii.kin[J], eii.occ[J])*widths[i];
    }
    tmp /= 1.4142; // Electron mass = 1 in atomic units
    return tmp;
}

double Distribution::tbr_int (const CustomDataType::EIIdata& eii, const int J) const {
    // Detailed-Balance dual of eii per atom for data corresponding to idx i.
    // TODO: check that this integral makes SENSE

    double tmp = 0;
    for (size_t i = 0; i < size; i++) {
        double x = 0;
        for (size_t j = 0; j < size; j++){
            double y= widths[j]*2*grid[j]*Dipole::DsigmaBEB( grid[i], grid[j], eii.ionB[J], eii.kin[J], eii.occ[J]);
            y *= f[i_from_e(grid[j] - grid[i] - eii.ionB[J])];
            x+= y;
        }
        tmp += x * f[i]*widths[i];
    }
    tmp *= 248.05021; // (2pi)^3
    return tmp;
}

void Distribution::set_maxwellian(double N, double T){
    T *= Constant::kb_in_au; // convert T to units of Ha
    for (size_t i=0; i<size; i++){
        f[i] = N*exp(-grid[i]/T)/pow(grid[i]*Constant::Pi*T, 0.5);
    }
}

std::string Distribution::get_energies(std::string units){
    std::stringstream ss;
    double conv=1;
    if (units == "eV") conv=Constant::eV_in_au;
    else if (units == "J") conv=Constant::eV_in_au*1.60217662e-19;
    else if (units == "Ha") conv=1;
    else if (units == "Rd") conv=2;

    for (auto& e : grid){
        ss<<e*conv<<' ';
    }
    return ss.str();
}

// Taken verbatim from Rockwood as quoted by Morgan and Penetrante in ELENDIF
void Distribution::add_Qee(const Distribution& F) {
    double alpha = 2.9619219588; // 2/3 pi e^4 * (2/m)^1/2, a.u.
    alpha *= 4.6; // Coulomb logarithm, estimated from ln(100) for Plasma parameter [very rough]
    double psi,dpsi,tmp; // psi = 3 I_0^e f(e') de' -1/e I_0^e e'f(e') de' + 2e^0.5 I_e^inf n(e') e'^-0.5 de'

    double integral_nde = 0;
    double integral_nede = 0;
    double integral_nsde = 0; // integral n(ep)/sqrt(ep)

    // Number of beginning and end points to exclude for numerical derivatives
    const int der_order = 1;
    double g[size]; // define g(e) = dn/de -n/2e

    // Nasty numerical differentiation
    for (size_t i=der_order; i<der_order*2; i++){
        double h1 = widths[i-1];
        double h2 = widths[i];
        double de1 = (F.f[i] - F.f[i-1]) /h1;
        double de2 = (F.f[i+1] - F.f[i]) /h2;
        g[i] = (h2*de1 + h1*de2) / (h1+h2) - 2 * F.f[i]/grid[i];

        // Compute edge values of f ignoring 2nd order gradient terms...
        this->f[i] = 3*F.f[i]*F.f[i];
        this->f[i] += psi*g[i];
        this->f[i] *= alpha/pow(grid[i],0.5);
    }

    for (size_t i=size-der_order*2; i<size-der_order; i++){
        double h1 = widths[i-1];
        double h2 = widths[i];
        double de1 = (F.f[i] - F.f[i-1]) /h1;
        double de2 = (F.f[i+1] - F.f[i]) /h2;
        g[i] = (h2*de1 + h1*de2) / (h1+h2) - 2 * F.f[i]/grid[i];

        // Compute edge values of f ignoring 2nd order gradient terms...
        this->f[i] = 3*F.f[i]*F.f[i];
        this->f[i] += psi*g[i];
        this->f[i] *= alpha/pow(grid[i],0.5);
    }

    for (size_t i=0; i < der_order; i++){
        // Compute edge values of f ignoring gradient terms...
        this->f[i] = alpha*3*F.f[i]*F.f[i]/pow(grid[i],0.5);
        size_t j=size-der_order+i;
        this->f[j] = alpha*3*F.f[j]*F.f[j]/pow(grid[j],0.5);
    }

    for (size_t i=der_order*2; i < size-der_order*2; i++){
        double e = grid[i];
        double e05 = pow(e,0.5);
        tmp = F.f[i]*widths[i];

        integral_nde += tmp;
        integral_nede += e*tmp;
        integral_nsde += tmp/e05;

        dpsi = (integral_nede/e/e + integral_nsde/e05);
        psi = 3*integral_nde -integral_nede/e + 2*e05*integral_nsde;

        double h1 = widths[i-1];
        double h2 = widths[i];
        double de1 = (g[i] - g[i-1]) /h1;
        double de2 = (g[i+1] - g[i]) /h2;
        tmp = (h2*de1 + h1*de2) / (h1+h2) - 2 * g[i]/grid[i];

        this->f[i] = 3*F.f[i]*F.f[i]/e05;
        this->f[i] += 2*e05*e05*e05*tmp;
        this->f[i] += psi*g[i]/e05;

        this->f[i] *= alpha;
    }
}



// f geometry:
/*
                 ___________________
        ________|                   |
   f[0]|  f[1]  |      f[2]         |
   w[0]|  w[1]  |      w[2]         |
  |    |        |                   |
  r[0] r[1]     r[2]                r[3]

*/

double Distribution::e_from_i(size_t i){
    // Linear grid
    return min_e + i*(max_e - min_e)/size;
}

size_t Distribution::i_from_e(double e){
    long i= std::floor((e - min_e) * size /(max_e - min_e));
    if (i < 0) return 0;
    if (i >= size) return size-i;
    return i;
}

void Distribution::addDeltaLike(double e, double mass){
    size_t i = i_from_e(e);
    f[i] += mass/widths[i];
}

state_type::state_type(){
    atomP.resize(P_sizes.size());
    for (size_t i = 0; i < atomP.size(); i++) {
        atomP[i].resize(P_sizes[i]);
    }
}


// Critical vector-space devices
state_type& state_type::operator+=(const state_type &s){
    for (size_t r = 0; r < atomP.size(); r++) {
        for (size_t i = 0; i < atomP[r].size(); i++) {
            atomP[r][i] += s.atomP[r][i];
        }
    }
    F += s.F;
    return *this;
}

state_type& state_type::operator*=(const double x){
    for (size_t r = 0; r < atomP.size(); r++) {
        for (size_t i = 0; i < atomP[r].size(); i++) {
            atomP[r][i] *= x;
        }
    }
    F *= x;
    return *this;
}

// state_type state_type::operator+(const state_type& s2){
//     state_type retval = *this;
//     retval += s2;
//     return retval;
// }
//
// state_type state_type::operator*(double x){
//     state_type retval = *this;
//     retval *= x;
//     return retval;
// }

// convenience members
state_type& state_type::operator=(const double x){
    for (auto& P : atomP){
        for (auto& p : P) {
            p=x;
        }
    }
    F = 0;
    return *this;
}

// Resizes the container to fit all of the states present in the atom ensemble
void state_type::set_P_shape(const vector<RateData::Atom>& atomsys) {
    P_sizes.resize(atomsys.size());
    // make the P's the right size lmao
    for (size_t a = 0; a < atomsys.size(); a++) {
        P_sizes[a] = atomsys[a].num_conf;
    }
}


// Intended usage: cout<<s.atomP[a]<<endl;
ostream& operator<<(ostream& os, const state_type::bound_t& bound){
    for (size_t i=0; i<bound.size(); i++){
        os << bound[i] << " ";
    }
    return os;
}

ostream& operator<<(ostream& os, const Distribution& dist){
    for (size_t i= 0; i < Distribution::size; i++) {
        os<<dist.f[i]<<" ";
    }
    return os;
}

ostream& operator<<(ostream& os, const state_type& st){
    for (size_t a=0; a<st.atomP.size(); a++){
        os << st.atomP[a];
        if (a != st.atomP.size()-1)
            os<<"| ";
    }
    os << st.F;
    return os;
}
