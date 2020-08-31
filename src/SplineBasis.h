#ifndef SPLINEBASIS_CXX_H
#define SPLINEBASIS_CXX_H

#include <vector>
#include <eigen3/Eigen/SparseLU>
// #include <eigen3/Eigen/LU>
#include <iostream>

class BasisSet{
public:
    BasisSet(){};
    void set_parameters(size_t n, double min, double max);
    Eigen::VectorXd Sinv(const Eigen::VectorXd& deltaf);
    double operator()(size_t i, double x) const;
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
    size_t num_funcs;
    const static int BSPLINE_ORDER =2; // Cubics should be fine, increase if needed
protected:
    // Eigen::PartialPivLU<Eigen::MatrixXd > linsolver;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >  linsolver;
    std::vector<double> knot;
    double overlap(size_t j, size_t k) const;
};


#endif /* end of include guard: SPLINEBASIS_CXX_H */
