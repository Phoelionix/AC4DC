#include "RateSystem.h"

vector<double> Distribution::grid   = vector<double>(0);
vector<double> Distribution::widths = vector<double>(0);
vector<size_t> state_type::P_sizes  = vector<size_t>(0);
size_t Distribution::size=0;

void Distribution::set_elec_points(size_t n, double min_e, double max_e){
    size=n;
    grid.resize(n);
    widths.resize(n);
    double step = (max_e - min_e) / n;
    grid[0]=min_e;
    for (size_t i=1; i<n; i++){
        grid[i] = min_e + step*i;
        widths[i] = grid[i]-grid[i-1];
    }
    // Questionable decision: "Practical infinity"
    widths[n-1] = 2*max_e;
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

state_type state_type::operator+(const state_type& s2){
    state_type retval = *this;
    retval += s2;
    return retval;
}

state_type state_type::operator*(double x){
    state_type retval = *this;
    retval *= x;
    return retval;
}

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
