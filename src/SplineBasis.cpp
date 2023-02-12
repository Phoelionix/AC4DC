/**
 * @file SpineBasis.cpp
 * @brief @copybrief SplineBasis.h
 * 
 */
/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================*/

#include "SplineBasis.h"
#include <assert.h>
#include "Constant.h"
#include "BSpline.h"
#include <algorithm>

// classic false position method
double find_root(std::function<double(double)> func, double a, double b, double tolerance=1e-8){
    assert(b > a);
    assert(tolerance > 0);
    double x = 0;
    while (fabs(func(x)) > tolerance){
        // x = (a*func(b) - b*func(a))/(func(b) - func(a));
        x = (a+b)/2;
        // split interval into [a, x]u[x, b]
        double fx = func(x);
        if (fx*func(a) < 0){
            // root in interval [x, a]
            b = x;
        } else if (fx*func(b) < 0){
            a = x;
        } else {
            throw runtime_error("Root Find Error: f does not change sign");
        }
    }
    return x;
}


//                  value to find, sorted array, length of array, index that *arr refers to
// Returns an index idx such that val lies in the interval [arr[idx], arr[idx+1])
int r_binsearch(double val, const double* arr, int len, int idx) {
    if (len==1) return idx;
    int mid = len/2;
    // arr has indices i i+1 ... i+len-1
    if (*(arr + mid) > val) {
        return r_binsearch(val, arr, mid, idx);
    } else {
        return r_binsearch(val, arr+mid, len - mid, idx+mid);
    }
}

/// @todo make variable names clearer
int lower_idx(double val, const std::vector<double>& arr){
    return r_binsearch(val, arr.data(), arr.size(), 0);
}

/// @todo make variable names clearer
int BasisSet::i_from_e(double e) {
    // size_t idx = r_binsearch(e, avg_e.data(), avg_e.size(), 0);
    size_t idx = lower_idx(e, avg_e);
    if (idx == avg_e.size()) return 0;
    // returns the closer one
    return (avg_e[idx+1] - e > e - avg_e[idx]) ? idx : idx +1;
}

int BasisSet::lower_i_from_e(double e) {
    // size_t idx = r_binsearch(e, avg_e.data(), avg_e.size(), 0);
    return lower_idx(e, avg_e);
}

/**
 * @brief Constructs the grid over which the electron distribution is solved. 
 * @details Idea: Want num_funcs usable B-splines
 *  If grid has num_funcs+k+1 points, num_funcs of these are usable splines
 *  grid layout for open boundary:
 *  t0=t1=t2=...=tk-1,   tn+1 = tn+2 =... =tn+k
 *  homogeneous Dirichlet boundary:
 *  t0=t1=t2=...=tk-2,   tn+2 =... =tn+k
 *  Neumann boundary:
 *  t0=t1=...=tk-3, tn+3=...=tn+k
 *     
 *  boundary at minimm energy enforces energy conservation
 * 
 * @param gt grid type, also referred to as grid_style.
 */
void BasisSet::set_knot(const GridSpacing& gt){
    int Z_0 = gt.zero_degree_0;
    int Z_inf = gt.zero_degree_inf;
    if( BSPLINE_ORDER - Z_0 < 0){
        std::cerr <<"Failed to set boundary condition at 0, defaulting to "<<0<<std::endl;
        Z_0 = 0;
    }
    if (BSPLINE_ORDER - Z_inf < 0){
        std::cerr <<"Failed to set boundary condition at infinity, defaulting to "<<BSPLINE_ORDER<<std::endl;
        Z_inf = BSPLINE_ORDER;
    } 

    size_t start = BSPLINE_ORDER-Z_0-1;
    size_t num_int = num_funcs + Z_inf - start;
    
 
    double A_lin = (_max - _min)/(num_int-1);
    double A_sqrt = (_max - _min)/(num_int-1)/(num_int-1);
    // exponential grid
    if ( (gt.mode == GridSpacing::exponential) && _min <= 0) {
        throw runtime_error("Cannot construct an exponential grid with zero minimum energy.");
    }
    double A_exp = _min;
    double lambda_exp = (log(_max) - log(_min))/(num_int-1);
    // hybrid linear-linear grid

    size_t M = gt.num_low + 1;

    assert(num_funcs > M); // TODO should only give a warning if not relevant to grid type. -S.P.
    assert(_max > gt.transition_e);

    double B_hyb = (gt.transition_e - _min)/M;
    double C_hyb = (_max - gt.transition_e)/(num_funcs - M);
    
    double p_powlaw = (log(_max-_min) - log(gt.transition_e - _min))/(log(1.*num_funcs/M));
    double A_powlaw = (_max - _min)/pow(num_funcs, p_powlaw);
    
    double p_powlaw_test = (log(_max-_min) - log(2.51 - _min))/(log(1.*num_funcs/M));   // transition e is in Ha, which are ~27 * larger than eV.
    double A_powlaw_test = (_max - _min)/pow(num_funcs, p_powlaw_test);

    /// Wiggle destroyer.  Continuous piecewise.
    // e.g. 6k eV photon beam. Have delta-like stuff in MB and dirac.
    // Idea: Define four energy regions: (eV) <-0--low--10--mid--200-trans.--2000--high--10000--surplus->
    // We want to have higher density in the mid and high energy regions, and low in the low, transition, and surplus regions.
    // Start with power law since, numerical analysis aside, it nearly works already.
    // 

    std::vector<double> hyb_powlaw_power (_region_bndry_index.size() - 1,0.);
    std::vector<double> hyb_powlaw_factor (_region_bndry_index.size() - 1,0.);
    // Params that define region boundaries:
    int n_M, n_N; //Index of first point in region/next region.  
    double E_M, E_N; // Corresponding energies.
    // R region boundaries and 1 power law for each region --> R - 1 power laws.
    for (size_t i=0; i< hyb_powlaw_power.size(); i++){
        n_M = _region_bndry_index[i];
        E_M = _region_bndry_energy[i];
        n_N = _region_bndry_index[i+1];
        E_N = _region_bndry_energy[i+1];
        
        hyb_powlaw_power[i] = (log(E_N - _min) - log(E_M - _min))/(log(1.*n_N) - log(1.*n_M));
        hyb_powlaw_factor[i] = (E_N - _min)/pow(n_N, hyb_powlaw_power[i]);
    }

    std::cout<<"[ Knot ] Using splines of "<<BSPLINE_ORDER<<"th order "<<std::endl;
    std::cout<<"[ Knot ] (i.e. piecewise polynomial has leading order x^"<<BSPLINE_ORDER-1<<")"<<std::endl;
    if ( gt.mode == GridSpacing::powerlaw) {
        std::cout<<"[ Knot ] Using power-law exponent "<<p_powlaw<<std::endl;
    }
    if ( gt.mode == GridSpacing::hybrid_power){
        std::cout<<"[ Knot ] Using power-law exponents: "; 
        for (int j = 0; j < hyb_powlaw_power.size(); j++){
            cout << hyb_powlaw_power[j] << ", ";
        }
        std::cout << std::endl; 
        std::cout<<"[ Knot ] Using power-law factors: ";
        for (int j = 0; j < hyb_powlaw_factor.size(); j++){
            cout << hyb_powlaw_factor[j] << ", ";
        }
        std::cout << std::endl;         
    }


    // std::function<double(double)> f = [=](double B) {return _min*exp(B*M/(_max - B*(num_funcs-M))) + B*(num_funcs-M) - _max;};
    // double B_hyb = 0;
    // if (gt.mode == GridSpacing::hybrid){
    //     B_hyb = find_root(f, 0, _max/(num_funcs-M + 1));
    // }
    // double C_hyb = _max - B_hyb * num_funcs;
    // double lambda_hyb = (log(B_hyb*M+C_hyb) - log(_min))/M;

    knot.resize(num_funcs+BSPLINE_ORDER+1); // num_funcs == n + 1

    // Set the k - zero_deg repeated knots
    
    for (size_t i=0; i<start; i++) {
        knot[i] = _min;
    }

    // TODO:
    // - Move this grid-construction code to GridSpacing.hpp
    // - Refactor to incorporate boundaries max and min into the GridSpacing object
    // - Unfuck the boundary conditions (seems to break for BSPLINE_ORDER=3?)

    // Keep in mind: At i=0, the value of knot will still be _min
    for(size_t i=start; i<=num_funcs + Z_inf; i++) {
        switch (gt.mode)
        {
        case GridSpacing::linear:
            knot[i] = _min + A_lin*(i-start); // linear spacing
            break;
        case GridSpacing::quadratic:
            knot[i] = _min + A_sqrt*(i-start)*(i-start); // square root spacing
            break;
        case GridSpacing::exponential:
            knot[i] = A_exp*exp(lambda_exp*(i-start));
            break;
        case GridSpacing::hybrid:
            // knot[i] = (i >=M ) ? B_hyb*(i-start) + C_hyb : A_exp * exp(lambda_hyb*(i-start));
            knot[i] = (i-start < M) ? B_hyb*(i-start) + _min : C_hyb*(i-start-M) + gt.transition_e;
            break;
        case GridSpacing::powerlaw:
            knot[i] = A_powlaw * pow(i-start, p_powlaw) + _min;
            cout << knot[i] * Constant::eV_per_Ha  << "," << i << endl;
            break;
        case GridSpacing::test:
            knot[i] = A_powlaw_test * pow(i - start, p_powlaw_test) + _min;
            break;
        case GridSpacing::hybrid_power:
        {
            //  Get the index of the region that this point is part of.
            size_t rgn = 0;
            for( ; rgn < hyb_powlaw_power.size(); rgn++){
                if(rgn == hyb_powlaw_power.size() - 1
                 ||i - start < _region_bndry_index[rgn+1]){
                    break;
                 }
            }
            knot[i] = hyb_powlaw_factor[rgn] * pow(i - start, hyb_powlaw_power[rgn]) + _min;
            cout << knot[i] * Constant::eV_per_Ha  << "," << rgn << ", i: " << i << endl;
            break;
        }
        default:
            throw std::runtime_error("Grid spacing has not been defined.");
            break;
        }
        
    }

    // t_{n+1+z_infinity} has been set now. Repeat it until the end.
    for (size_t i= num_funcs + 1 + Z_inf; i< num_funcs+BSPLINE_ORDER; i++) {
        knot[i] = knot[num_funcs + Z_inf];
    }
    

    #ifdef DEBUG
    std::cerr<<"Knot: [ ";
    for (auto& ki : knot) {
        std::cerr<<ki<<" ";
    }
    std::cerr<<"]\n";
    #endif
}

/**
 * @brief Sets up the B-spline knot to have the appropriate shape (respecting boundary conditions)
 * @details 
 * 
 * @param _num_funcs the number of B-splines to use in the basis. 
 * @param min minimum energy
 * @param max maximum energy
 * @param gt gt.zero_degree_0: The number of derivatives to set to zero: 0 = open conditions, 1=impose f(0)=0, 2=impose f(0)=f'(0) =0
 * @param  
 */
void BasisSet::set_parameters(size_t num_of_funcs, double min, double max, const GridSpacing& gt, std::vector<int> region_bndry_index, std::vector<double> region_bndry_energy) {
    this->_min = min;
    this->_max = max;
    this->_region_bndry_index = region_bndry_index;   // TODO: refactor, can replace min and max with array of region start and region ends. -S.P.
    this -> _region_bndry_energy = region_bndry_energy;   

    num_funcs = num_of_funcs;
    // Override params
    if(gt.mode == GridSpacing::hybrid_power){
        num_funcs = _region_bndry_index[_region_bndry_index.size()-1];
        this->_max = _region_bndry_energy[_region_bndry_energy.size()-1];
    }
    // Idea: Want num_funcs usable B-splines
    // If grid has num_funcs+k+1 points, num_funcs of these are usable splines (num_int = num_funcs + Z_inf - start)
    // grid layout for open boundary:
    // t0=t1=t2=...=tk-1,   tn+1 = tn+2 =... =tn+k
    // homogeneous Dirichlet boundary:
    // t0=t1=t2=...=tk-2,   tn+2 =... =tn+k
    // Neumann boundary:
    // t0=t1=...=tk-3, tn+3=...=tn+k
    
    // boundary at minimm energy enforces energy conservation
    set_knot(gt);
    

    
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(BSPLINE_ORDER*2*num_funcs);
    
    // Eigen::MatrixXd S(num_funcs, num_funcs);
    for (size_t i=0; i<num_funcs; i++) {
        for (size_t j=i+1; j<num_funcs; j++) {
            double tmp=overlap(i, j);
            if (tmp != 0) {
                tripletList.push_back(Eigen::Triplet<double>(i,j,tmp));
                tripletList.push_back(Eigen::Triplet<double>(j,i,tmp));
                // S(i,j) = tmp;
                // S(j,i) = tmp;
            }
        }
        // S.insert(i,i) = overlap(i,i);
        tripletList.push_back(Eigen::Triplet<double>(i,i,overlap(i, i)));
        // S(i,i) = overlap(i,i);
    }
    /// Compute overlap matrix   //s\/ If we use the grid basis (with num_func splines) the splines overlap. S takes us to a spline-independent basis I believe.
    Eigen::SparseMatrix<double> S(num_funcs, num_funcs);  //s\/ rows and cols. 
    S.setFromTriplets(tripletList.begin(), tripletList.end());
    // Precompute the LU decomposition
    linsolver.analyzePattern(S);
    linsolver.isSymmetric(true);
    linsolver.factorize(S);
    if(linsolver.info()!=Eigen::Success) {
        throw runtime_error("Factorisation of overlap matrix failed!");
    }

    avg_e.resize(num_funcs);
    log_avg_e.resize(num_funcs);
    areas.resize(num_funcs);
    for (size_t i = 0; i < num_funcs; i++) {
        // Chooses the 'center' of the B-spline
        avg_e[i] = (this->supp_max(i) + this->supp_min(i))/2 ;
        log_avg_e[i] = log(avg_e[i]);
        double diff = (this->supp_max(i) - this->supp_min(i))/2;
        // Widths stores the integral of the j^th B-spline
        areas[i] = 0;
        for (size_t j=0; j<10; j++){
            areas[i] += gaussW_10[j]*(*this)(i, gaussX_10[j]*diff+ avg_e[i]);
        }
        areas[i] *= diff;
    }
}

/**
 * @brief Changes basis of vector deltaf from spline basis to grid basis, see details.
 * @details From SplineBasis.cpp - "If grid has num_funcs+k+1 points, num_funcs of these are usable splines".
 */
Eigen::VectorXd BasisSet::Sinv(const Eigen::VectorXd& deltaf) {
    // Solves the linear system S fdot = deltaf
    return linsolver.solve(deltaf);
}

/**
 * @brief Changes basis of matrix J from spline basis to grid basis, see details.
 * @details From SplineBasis.cpp - "If grid has num_funcs +k+1 points, num_funcs of these are usable splines".
 */
Eigen::MatrixXd BasisSet::Sinv(const Eigen::MatrixXd& J) {
    // Solves the linear system S fdot = deltaf
    return linsolver.solve(J);
}

double BasisSet::raw_bspline(size_t i, double x) const {
    assert(i < num_funcs && i >= 0);
    return BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
}

double BasisSet::raw_Dbspline(size_t i, double x) const {
    assert(i < num_funcs && i >= 0);
    return BSpline::DBSpline<BSPLINE_ORDER>(x, &knot[i]);
}

double BasisSet::operator()(size_t i, double x) const {
    assert(i < num_funcs && i >= 0);
    // Returns the i^th B-spline of order BSPLINE_ORDER
    
    static_assert(USING_SQRTE_PREFACTOR);
    return pow(x,0.5)*BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
    // return BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
    
}

// Returns the i^th B-spline of order BSPLINE_ORDER evaluated at x (safe version)
double BasisSet::at(size_t i, double x) const {
    if(i >= num_funcs || i<0) return 0;
    return operator()(i, x);
}

// Returns the first derivative of the i^th B-spline of order BSPLINE_ORDER evaluated at x
double BasisSet::D(size_t i, double x) const {
    assert(i < num_funcs && i >= 0);
    static_assert(BSPLINE_ORDER > 1);
    
    static_assert(USING_SQRTE_PREFACTOR);
    return pow(x,0.5)*BSpline::DBSpline<BSPLINE_ORDER>(x, &knot[i]) + 0.5*pow(x,-0.5)*BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
    
    // return BSpline::DBSpline<BSPLINE_ORDER>(x, &knot[i]);
    
}

double BasisSet::overlap(size_t j,size_t i) const{
    // computes the overlap of the jth raw bspline with the ith modified bspline
    double b = min(supp_max(i), supp_max(j));
    double a = max(supp_min(i), supp_min(j));
    if (a >= b) return 0;
    double tmp=0;
    for (int m=0; m<10; m++) {
        double e = gaussX_10[m]*(b-a)/2 + (b+a)/2;
        tmp += gaussW_10[m]*raw_bspline(j, e)*(*this)(i, e);
        // tmp += gaussW_10[m]*(*this)(j, e)*(*this)(i, e);
    }
    tmp *= (b-a)/2;
    return tmp;
}

// returns a vector of intervals (bottom, K1), (K1, K2), (K2, top)
// such that the splines are polynomial on these intervals
std::vector<std::pair<double,double>> BasisSet::knots_between(double bottom, double top) const {
    std::vector<std::pair<double,double> > intervals(0);
    // NOTE:
    // These calls can and should should be replaced with std::lower_bound calls.
    if (bottom >= top) return intervals;
    size_t bottomidx = lower_idx(bottom, knot);
    size_t topidx = lower_idx(top, knot) + 1;
    if (topidx >= knot.size()) topidx = knot.size() -1;
    for (size_t i = bottomidx; i<= topidx; i++){
        intervals.push_back(std::pair<double,double>(knot[i],knot[i+1]));
    }
    return intervals;
}