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
    assert(zero_degree < BSPLINE_ORDER -1 && "Failed to set boundary consition");

    for (int i=0; i<BSPLINE_ORDER-zero_degree; i++) {
        knot[i] = _min;
    }

    // all derivatives vanish at "infinity": no special treatment
    int n = num_funcs;
    

    double A_lin = (_max - _min)/(n-1);
    double A_sqrt = (_max - _min)/(n-1)/(n-1);
    // exponential grid
    if ( (gt == GridSpacing::exponential || gt == GridSpacing::hybrid) && _min <= 0) {
        throw runtime_error("Cannot construct an exponential grid with zero minimum energy.");
    }
    double A_exp = _min;
    double lambda_exp = (log(_max) - log(_min))/(n-1);
    // hybrid exponentio-linear grid
    int M = n/10;
    double max_step = A_lin * 10;
    std::function<double(double)> f = [=](double l) {return exp(l*M) - exp(l*(M-1)) - max_step/_min; };
    double lambda_hyb = find_root(f, 0, lambda_exp*100);
    double B_hyb = (A_exp*exp(lambda_hyb*M) - _max)/(M-n);
    double C_hyb = _max - B_hyb * n;

    for(int i=BSPLINE_ORDER-zero_degree; i<n+BSPLINE_ORDER; i++) {
        switch (gt)
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
}

// Sets up the B-spline knot to have the appropriate shape (respecting boundary conditions)
void BasisSet::set_parameters(size_t n, double min, double max, int zero_degree, GridSpacing gt) {
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
        std::cerr<<"Factorisation of overlap matrix failed!"<<endl;
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
