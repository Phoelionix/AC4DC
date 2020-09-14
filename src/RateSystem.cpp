#include "RateSystem.h"
#include "Dipole.h"
#include <math.h>
// #include <stringstream>
// #define NDEBUG


vector<size_t> state_type::P_sizes  = vector<size_t>(0);


state_type::state_type() {
    atomP.resize(P_sizes.size());
    for (size_t i = 0; i < atomP.size(); i++) {
        atomP[i].resize(P_sizes[i]);
    }
}


// Critical vector-space operators
state_type& state_type::operator+=(const state_type &s) {
    for (size_t r = 0; r < atomP.size(); r++) {
        for (size_t i = 0; i < atomP[r].size(); i++) {
            atomP[r][i] += s.atomP[r][i];
        }
    }
    F += s.F;
    return *this;
}

state_type& state_type::operator*=(const double x) {
    for (size_t r = 0; r < atomP.size(); r++) {
        for (size_t i = 0; i < atomP[r].size(); i++) {
            atomP[r][i] *= x;
        }
    }
    F *= x;
    return *this;
}

// convenience members
state_type& state_type::operator=(const double x) {
    for (auto& P : atomP) {
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

// Returns the L1 norm
double state_type::norm() const {
    double n = 0;
    for (auto& P : atomP) {
        for (auto& p : P) {
            n += fabs(p);
        }
    }
    n += F.norm();
    return n;
}


// Intended usage: cout<<s.atomP[a]<<endl;
ostream& operator<<(ostream& os, const bound_t& bound) {
    for (size_t i=0; i<bound.size(); i++) {
        os << bound[i] << " ";
    }
    return os;
}

ostream& operator<<(ostream& os, const state_type& st) {
    for (size_t a=0; a<st.atomP.size(); a++) {
        os << st.atomP[a];
        if (a != st.atomP.size()-1)
            os<<"| ";
    }
    os << st.F;
    return os;
}
