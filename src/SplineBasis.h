#ifndef SPLINEBASIS_CXX_H
#define SPLINEBASIS_CXX_H

#include <vector>
#include <eigen3/Eigen/SparseLU>
// #include <eigen3/Eigen/LU>
#include <iostream>
#include "GridSpacing.hpp"


class BasisSet{
public:
    BasisSet() {};
    void set_parameters(size_t nfuncs, double min, double max, GridSpacing gt, int zero_degree);
    Eigen::VectorXd Sinv(const Eigen::VectorXd& deltaf);
    // Returns the ith basis function evaluated at point x
    double operator()(size_t i, double x) const;
    // Returns the first derivative of the ith basis function at point x
    double D(size_t i, double x) const;
    double at(size_t i, double x) const;
    inline double supp_max(unsigned i) const{
        return knot[i+BSPLINE_ORDER];
    }
    inline double supp_min(unsigned i) const{
        return knot[i];
    }
    size_t gridlen() const{
        return knot.size();
    }
    double grid(size_t i) const{
        return knot[i];
    }
    double min_elec_e() {return _min;};
    double max_elec_e() {return _max;};
    size_t num_funcs;
    const static int BSPLINE_ORDER = 3; // 1 = rectangles, 2=linear, 3=quadratic
    std::vector<double> avg_e;
    std::vector<double> avg_e_sqrt;
    std::vector<double> widths;
protected:
    // Eigen::PartialPivLU<Eigen::MatrixXd > linsolver;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >  linsolver;
    std::vector<double> knot;
    double overlap(size_t j, size_t k) const;
    void set_knot(unsigned zero_degree, GridSpacing gt);
    

    double _min;
    double _max;
};


#endif /* end of include guard: SPLINEBASIS_CXX_H */
