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
void BasisSet::set_knot(int zero_degree, GridSpacing gt){
    assert( BSPLINE_ORDER - zero_degree >=0 && "Failed to set boundary condition");

    for (int i=0; i<BSPLINE_ORDER-zero_degree; i++) {
        knot[i] = _min;
    }

    // all derivatives vanish at "infinity": no special treatment
    int n = num_funcs;
    

    double A_lin = (_max - _min)/(n-1);
    double A_sqrt = (_max - _min)/(n-1)/(n-1);
    // exponential grid
    if ( (gt.mode == GridSpacing::exponential || gt.mode == GridSpacing::hybrid) && _min <= 0) {
        throw runtime_error("Cannot construct an exponential grid with zero minimum energy.");
    }
    double A_exp = _min;
    double lambda_exp = (log(_max) - log(_min))/(n-1);
    // hybrid exponentio-linear grid

    
    int M = gt.num_exp;
    
    
    std::function<double(double)> f = [=](double B) {return _min*exp(B*M/(_max - B*(n-M))) + B*(n-M) - _max;};
    double B_hyb = find_root(f, 0, _max/(n-M + 1));
    
    double C_hyb = _max - B_hyb * n;
    double lambda_hyb = (log(B_hyb*M+C_hyb) - log(_min))/M;
    
    
    for(int i=BSPLINE_ORDER-zero_degree; i<n+BSPLINE_ORDER; i++) {
        switch (gt.mode)
        {
        case GridSpacing::linear:
            knot[i] = _min + A_lin*i; // linear spacing
            break;
        case GridSpacing::quadratic:
            knot[i] = _min + A_sqrt*i*i; // square root spacing
            break;
        case GridSpacing::exponential:
            knot[i] = A_exp*exp(lambda_exp*i);
            break;
        case GridSpacing::hybrid:
            knot[i] = (i >=M ) ? B_hyb*i + C_hyb : A_exp * exp(lambda_hyb*i);
            break;
        default:
            throw std::runtime_error("Grid spacing has not been defined.");
            break;
        }
        
    }

    #ifdef DEBUG
    std::cerr<<"Knot: [";
    for(int i=0; i<n+BSPLINE_ORDER; i++) {
        std::cerr<<knot[i]<<" ";
    }
    std::cerr<<"]\n";
    #endif
}

// Sets up the B-spline knot to have the appropriate shape (respecting boundary conditions)
void BasisSet::set_parameters(size_t n, double min, double max, GridSpacing gt, int zero_degree) {
    // zero_degree: The number of derivatives to set to zero: 0 = open conditions, 1=impose f(0)=0, 2=impose f(0)=f'(0) =0
    // n: the number of Bsplins to use in the basis
    // min: minimum energy
    // max: maximum energy
    this->_min = min;
    this->_max = max;

    num_funcs = n;
    // Idea: Want n usable B-splines
    // If grid has n+k+1 points, n of these are usable splines
    // grid layout for open boundary:
    // t0=t1=t2=...=tk-1,   tn+1 = tn+2 =... =tn+k
    // homogeneous Dirichlet boundary:
    // t0=t1=t2=...=tk-2,   tn+2 =... =tn+k
    // Neumann boundary:
    // t0=t1=...=tk-3, tn+3=...=tn+k
    knot.resize(n+BSPLINE_ORDER);
    // boundary at minimm energy enforces
    set_knot(zero_degree, gt);

    avg_e.resize(n);
    avg_e_sqrt.resize(n);
    widths.resize(n);
    for (size_t i = 0; i < n; i++) {
        // Chooses the 'center' of the B-spline
        if (BSPLINE_ORDER%2==0){
            avg_e[i] = this->knot[i + BSPLINE_ORDER/2];
        } else  {
            avg_e[i] = (this->knot[i + BSPLINE_ORDER/2] + this->knot[i + BSPLINE_ORDER/2 + 1]) /2;
        }
        avg_e_sqrt[i] = pow(avg_e[i],0.5);
        double diff = (this->supp_max(i) - this->supp_min(i))/2;
        widths[i] = 0;
        for (int j=0; j<10; j++){
            widths[i] += gaussW_10[j]*(*this)(i, gaussX_10[j]*diff+ avg_e[i]);
        }
        widths[i] *= diff;
    }
    
    // Compute overlap matrix
    Eigen::SparseMatrix<double> S(num_funcs, num_funcs);
    // Eigen::MatrixXd S(num_funcs, num_funcs);
    for (int i=0; i<num_funcs; i++) {
        for (int j=i+1; j<num_funcs; j++) {
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
}

Eigen::VectorXd BasisSet::Sinv(const Eigen::VectorXd& deltaf) {
    // Solves the linear system S fdot = deltaf
    return linsolver.solve(deltaf);
}

double BasisSet::operator()(size_t i, double x) const{
    assert(i < num_funcs && i >= 0);
    // Returns the i^th B-spline of order BSPLINE_ORDER
    return BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
}

double BasisSet::at(size_t i, double x) const{
    if(i >= num_funcs || i<0) return 0;
    // Returns the i^th B-spline of order BSPLINE_ORDER
    return BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
}

double BasisSet::overlap(size_t j,size_t i) const{
    int diff = (int)j -  (int)i;
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
