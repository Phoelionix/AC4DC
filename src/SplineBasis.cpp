#include "SplineBasis.h"
#include <assert.h>
#include "Constant.h"
#include "BSpline.h"

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

// Constructs the grid over which the electron distribution is solved. 

void BasisSet::set_knot(const GridSpacing& gt){
    int Z_0 = gt.zero_degree_0;
    int Z_inf = gt.zero_degree_inf;
    if( BSPLINE_ORDER - Z_0 < 0){
        std::cerr <<"Failed to set boundary condition at 0, defaulting to "<<BSPLINE_ORDER<<std::endl;
        Z_0 = BSPLINE_ORDER;
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
    if ( (gt.mode == GridSpacing::exponential || gt.mode == GridSpacing::hybrid) && _min <= 0) {
        throw runtime_error("Cannot construct an exponential grid with zero minimum energy.");
    }
    double A_exp = _min;
    double lambda_exp = (log(_max) - log(_min))/(num_int-1);
    // hybrid exponentio-linear grid

    
    size_t M = gt.num_exp;
    
    
    std::function<double(double)> f = [=](double B) {return _min*exp(B*M/(_max - B*(num_int-M))) + B*(num_int-M) - _max;};
    
    double B_hyb = 0;
    if (gt.mode == GridSpacing::hybrid){
        B_hyb = find_root(f, 0, _max/(num_int-M + 1));
    }
    
    double C_hyb = _max - B_hyb * num_int;
    double lambda_hyb = (log(B_hyb*M+C_hyb) - log(_min))/M;
    

    knot.resize(num_funcs+BSPLINE_ORDER); // num_funcs == n + 1

    // Set the k - zero_deg repeated knots
    
    for (size_t i=0; i<start; i++) {
        knot[i] = _min;
    }

    // At i=0, the value of knot will still be _min
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
            knot[i] = (i >=M ) ? B_hyb*(i-start) + C_hyb : A_exp * exp(lambda_hyb*(i-start));
            break;
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

// Sets up the B-spline knot to have the appropriate shape (respecting boundary conditions)
void BasisSet::set_parameters(size_t num_int, double min, double max, const GridSpacing& gt) {
    // gt.zero_degree_0: The number of derivatives to set to zero: 0 = open conditions, 1=impose f(0)=0, 2=impose f(0)=f'(0) =0
    // num_int: the number of Bsplins to use in the basis
    // min: minimum energy
    // max: maximum energy
    this->_min = min;
    this->_max = max;

    num_funcs = num_int;
    // Idea: Want num_int usable B-splines
    // If grid has num_int+k+1 points, num_int of these are usable splines
    // grid layout for open boundary:
    // t0=t1=t2=...=tk-1,   tn+1 = tn+2 =... =tn+k
    // homogeneous Dirichlet boundary:
    // t0=t1=t2=...=tk-2,   tn+2 =... =tn+k
    // Neumann boundary:
    // t0=t1=...=tk-3, tn+3=...=tn+k
    
    // boundary at minimm energy enforces energy conservation
    set_knot(gt);
    
    // Compute overlap matrix
    Eigen::SparseMatrix<double> S(num_funcs, num_funcs);
    // Eigen::MatrixXd S(num_funcs, num_funcs);
    for (size_t i=0; i<num_funcs; i++) {
        for (size_t j=i+1; j<num_funcs; j++) {
            double tmp=overlap(i, j);
            if (tmp != 0) {
                S.insert(i,j) = tmp;
                S.insert(j,i) = tmp;
                // S(i,j) = tmp;
                // S(j,i) = tmp;
            }
        }
        S.insert(i,i) = overlap(i,i);
        // S(i,i) = overlap(i,i);
    }
    // Precompute the LU decomposition
    linsolver.analyzePattern(S);
    linsolver.isSymmetric(true);
    linsolver.factorize(S);
    if(linsolver.info()!=Eigen::Success) {
        throw runtime_error("Factorisation of overlap matrix failed!");
    }

    avg_e.resize(num_int);
    avg_e_sqrt.resize(num_int);
    areas.resize(num_int);
    for (size_t i = 0; i < num_int; i++) {
        // Chooses the 'center' of the B-spline
        avg_e[i] = (this->supp_max(i) + this->supp_min(i))/2 ;
        avg_e_sqrt[i] = pow(avg_e[i],0.5);
        double diff = (this->supp_max(i) - this->supp_min(i))/2;
        // Widths stores the integral of the j^th B-spline
        areas[i] = 0;
        for (size_t j=0; j<10; j++){
            areas[i] += gaussW_10[j]*(*this)(i, gaussX_10[j]*diff+ avg_e[i]);
        }
        areas[i] *= diff;
    }
}

Eigen::VectorXd BasisSet::Sinv(const Eigen::VectorXd& deltaf) {
    // Solves the linear system S fdot = deltaf
    return linsolver.solve(deltaf);
}

double BasisSet::operator()(size_t i, double x) const{
    assert(i < num_funcs && i >= 0);
    // Returns the i^th B-spline of order BSPLINE_ORDER
    
    // return (i==0) ? pow(x,0.5)*BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]) : BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
    // return pow(x,0.5)*BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
    return BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
    
}

// Returns the i^th B-spline of order BSPLINE_ORDER evaluated at x (safe version)
double BasisSet::at(size_t i, double x) const {
    if(i >= num_funcs || i<0) return 0;
    return operator()(i, x);
}

// Returns the dirst derivative of the i^th B-spline of order BSPLINE_ORDER evaluated at x
double BasisSet::D(size_t i, double x) const {
    assert(i < num_funcs && i >= 0);
    static_assert(BSPLINE_ORDER > 1);
    // if (i == 0){
    // return pow(x,0.5)*BSpline::DBSpline<BSPLINE_ORDER>(x, &knot[i]) + 0.5*pow(x,-0.5)*BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
    // } else {
    return BSpline::DBSpline<BSPLINE_ORDER>(x, &knot[i]);
    // }
    
}

double BasisSet::overlap(size_t j,size_t i) const{
    double b = min(supp_max(i), supp_max(j));
    double a = max(supp_min(i), supp_min(j));
    if (a >= b) return 0;
    double tmp=0;
    for (int m=0; m<10; m++) {
        double e = gaussX_10[m]*(b-a)/2 + (b+a)/2;
        tmp += gaussW_10[m]*(*this)(j, e)*(*this)(i, e);
    }
    tmp *= (b-a)/2;
    return tmp;
}
