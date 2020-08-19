#include "SplineBasis.h"
#include <assert.h>
#include "Constant.h"


namespace BSpline{
    // Template voodoo stolen from stackexchange
    template <unsigned k>
    double BSpline(double x, double *t)
    {
        if (*t <= x && x < *(t+k))
        {
            double a=0;
            double h = *(t+k-1)-*t;
            a += (fabs(h) > 1e-16) ? BSpline<k-1>(x, t) * (x - *t) / h : 0;
            h = (*(t+k) - *(t+1));
            a += (fabs(h) > 1e-16) ? BSpline<k-1>(x, (t+1)) * (*(t+k) - x) / h : 0;
            return a;
        }
        else
            return 0;
    }

    template <>
    double BSpline<1>(double x, double *t)
    {
        if (*t <= x && x < *(t+1))
            return 1.;
        else
            return 0.;
    }

};

// Sets up the B-spline knot to have the appropriate shape (respecting boundary conditions)
void BasisSet::set_parameters(size_t n, double min, double max) {
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
    // open boundary at minimm energy
    for (int i=0; i<BSPLINE_ORDER; i++){
        knot[i] = min;
    }
    double de = 1.*(max-min)/n;
    // all derivatives vanish at "infinity": no special treatment
    for(int i=BSPLINE_ORDER; i<n+BSPLINE_ORDER; i++){
        knot[i] = knot[i-1] + de;
    }
    // Compute overlap matrix
    Eigen::SparseMatrix<double> S(num_funcs, num_funcs);
    // Eigen::MatrixXd S(num_funcs, num_funcs);
    for (int i=0; i<num_funcs; i++){
        for (int j=i+1; j<num_funcs; j++){
            double tmp=overlap(i, j);
            if (tmp != 0){
                S.insert(i,j) = tmp;
                S.insert(j,i) = tmp;
                // S(i,j) = tmp;
                // S(j,i) = tmp;
            }
        }
        S.insert(i,i) = overlap(i,i);
        // S(i,i) = overlap(i,i);
    }
    // Compute the decomposition
    // linsolver.compute(S);
    linsolver.analyzePattern(S);
    linsolver.isSymmetric(true);
    // Compute the numerical factorization
    linsolver.factorize(S);
    if(linsolver.info()!=Eigen::Success) {
        std::cerr<<"Factorisation of overlap matrix failed!"<<endl;
    }
}

Eigen::VectorXd BasisSet::Sinv(const Eigen::VectorXd& deltaf){
    // Solves the linear system S fdot = deltaf
    return linsolver.solve(deltaf);
}

double BasisSet::operator()(size_t i, double x){
    assert(i < num_funcs);
    // Returns the i^th B-spline of order BSPLINE_ORDER
    return BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
}


double BasisSet::overlap(size_t j,size_t i){
    if (j<i){
        int s = i;
        i=j;
        j=s;
    }// j >= i
    if (j-i>=BSPLINE_ORDER) return 0;
    double b = supp_max(j);
    double a = supp_min(i);
    double tmp=0;
    for (int m=0; m<10; m++){
        double e = gaussX_10[m]*(b-a)/2 + (b+a)/2;
        tmp += gaussW_10[m]*(*this)(j, e)*(*this)(i, e);
    }
    tmp *= (b-a)/2;
    return tmp;
}
