#include "RateSystem.h"

vector<double> state_type::f_grid   = vector<double>(0);
vector<double> state_type::f_widths = vector<double>(0);
vector<size_t> state_type::P_sizes  = vector<size_t>(0);

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
    for (size_t i= 0; i < st.f.size(); i++) {
        os<<st.f[i]<<" ";
    }

    return os;
}

void state_type::set_elec_points(size_t n, double min_e, double max_e){
    f_grid.resize(n);
    f_widths.resize(n);
    double step = (max_e - min_e) / n;
    f_grid[0]=min_e;
    for (size_t i=1; i<n; i++){
        f_grid[i] = min_e + step*i;
        f_widths[i] = f_grid[i]-f_grid[i-1];
    }
    // Questionable decision: "Practical infinity"
    f_widths[n-1] = 2*max_e;
}

void state_type::print_info(){
    cout<<"P sizes:";
    for (auto& s : P_sizes){
        cout<<" "<<s;
    }
    cout<<endl<<"f description: "<<endl<<"e |";
    for (auto& s : f_grid) {
        cout<<" "<<s;
    }
    cout<<endl<<" w |";
    for (auto& s : f_widths) {
        cout<<" "<<s;
    }
    cout<<endl;
}
