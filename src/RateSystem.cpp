#include "RateSystem.h"
#include "Dipole.h"
#include <math.h>

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
            x += widths[j]*2*grid[j]*Dipole::DsigmaBEB( grid[i], grid[j], eii.ionB[J], eii.kin[J], eii.occ[J]) * f[grid[j] - grid[i] - eii.ionB[J]];
        }
        tmp += x * f[i]*widths[i];
    }
    tmp *= 248.05021; // (2pi)^3
    return tmp;
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
    for (auto& fi : F.f){
        fi=x;
    }
    return *this;
}

void state_type::print_info(){
    cout<<"P sizes:";
    for (auto& s : P_sizes){
        cout<<" "<<s;
    }
    cout<<endl<<"f description: "<<endl<<"e |";
    for (auto& s : Distribution::grid) {
        cout<<" "<<s;
    }
    cout<<endl<<" w |";
    for (auto& s : Distribution::widths) {
        cout<<" "<<s;
    }
    cout<<endl;
}

// Resizes the container to fit all of the states present in the atom ensemble
void state_type::set_P_shape(const vector<RateData::Atom>& atomsys) {
    P_sizes.resize(atomsys.size());
    // make the P's the right size lmao
    for (size_t a = 0; a < atomsys.size(); a++) {
        P_sizes[a] = atomsys[a].num_conf;
    }
}

ostream& operator<<(ostream& os, const state_type& st){
    os<<"[ ";
    for (size_t a=0; a<st.atomP.size(); a++){
        for (size_t i=0; i<st.atomP[a].size(); i++){
            os << st.atomP[a][i] << " ";
        }
        if (a != st.atomP.size()-1)
            os<<"| ";
    }
    os<<"] ";
    os<<st.F;

    return os;
}

ostream& operator<<(ostream& os, const Distribution& dist){
    for (size_t i= 0; i < Distribution::size; i++) {
        os<<dist.f[i]<<" ";
    }

    return os;
}
