#ifndef SPLINEBASIS_CXX_H
#define SPLINEBASIS_CXX_H

#include <vector>
#include <eigen3/Eigen/SparseLU>
// #include <eigen3/Eigen/LU>
#include <iostream>

enum class GridSpacing {linear, quadratic, exponential};

namespace{
    std::ostream& operator<<(std::ostream& os, GridSpacing gs) {
        switch (gs)
        {
        case GridSpacing::linear:
            os << "linear";
            break;
        case GridSpacing::quadratic:
            os << "quadratic";
            break;
        case GridSpacing::exponential:
            os << "exponential";
            break;
        default:
            os << "Unknown grid type";
            break;
        }
        return os;
    }

    std::istream& operator>>(std::istream& is, GridSpacing& gs) {
        std::string tmp;
        is >> tmp;
        if (tmp.length() == 0) {
            std::cerr<<"No grid type provided, defaulting to linear..."<<std::endl;
            gs = GridSpacing::linear;
            return is;
        }
        switch ((char) tmp[0])
        {
        case 'l':
            gs = GridSpacing::linear;
            break;
        case 'q':
            gs = GridSpacing::quadratic;
            break;
        case 'e':
            gs = GridSpacing::exponential;
            break;
        default:
            std::cerr<<"Unrecognised grid type \""<<tmp<<"\""<<std::endl;
            gs = GridSpacing::linear;
            break;
        }
        return is;
    }
}

class BasisSet{
public:
    BasisSet() {};
    void set_parameters(size_t nfuncs, double min, double max, int zero_degree, GridSpacing gt);
    Eigen::VectorXd Sinv(const Eigen::VectorXd& deltaf);
    double operator()(size_t i, double x) const;
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
protected:
    // Eigen::PartialPivLU<Eigen::MatrixXd > linsolver;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >  linsolver;
    std::vector<double> knot;
    double overlap(size_t j, size_t k) const;
    

    double _min;
    double _max;
};


#endif /* end of include guard: SPLINEBASIS_CXX_H */
