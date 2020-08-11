#ifndef RATESYSTEM_CXX_H
#define RATESYSTEM_CXX_H


#include <sstream>
#include <assert.h>
#include <iostream>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseCholesky>

#define BSPLINE_ORDER 4 // Cubics should be fine, increase if needed
#define GAUSS_ORDER 10

class BasisSet{
public:
    BasisSet(size_t n, double min, double max);
    Eigen::VectorXd Sinv(Eigen::VectorXd deltaf);
    double operator()(size_t i, double x);
    inline double supp_max(unsigned k);
    inline double supp_min(unsigned k);
    size_t num_funcs;
private:
    Eigen::SimplicalLLT<SparseMatrix<double> > cholmachine;
    std::vector<double> knot;
    double overlap(size_t j, size_t k);
    // Template voodoo stolen from stackexchange
    template <unsigned k>
    double BSpline(double x, double *t, unsigned n)
    {
        if (*t <= x && x < *(t+k))
        {
            double a = (x - *t) / (*(t+k-1) - *t);
            double b = (*(t+k) - x) / (*(t+k) - *(t+1));

            return a * BSpline<k-1>(x, t) + b * BSpline<k-1>(x, (t+1));
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

// Represents a statistical distribution of electrons. Internal units are atomic units.
class Distribution
{
public:
    Distribution() {
        f.resize(size);
    }

    double operator[](size_t n){
        return this->f[n];
    }

    double operator()(double e){
        return this->f[i_from_e(e)];
    }

    // vector-space algebra
    Distribution& operator+=(const Distribution& d){
        for (size_t i=0; i<size; i++) {
            f[i] += d.f[i];
        }
        return *this;
    }

    Distribution& operator*=(double x){
        for (auto& fi : f){
            fi *= x;
        }
        return *this;
    }

    Distribution& operator=(const Distribution& d){
        f = d.f;
        return *this;
    }

    Distribution& operator=(double y){
        for (auto& fi : f){
            fi=y;
        }
        return *this;
    }

    template<typename T>
    T expect(T (& g) (double e)){
        T tmp = 0;
        for (size_t i = 0; i < size; i++) {
            tmp += g(grid[i])*f[i]*widths[i];
        }
        return tmp;
    }

    static void Gamma_eii(Eigen::SparseMatrix<double>& Gamma, const CustomDataType::EIIdata& eii, size_t K, int k);
    static void Gamma_tbr(Eigen::SparseMatrix<double>& Gamma, const CustomDataType::EIIdata& eii, size_t J, size_t K, int k);
    void add_Qeii (size_t a, const Distribution& F, const bound_t& P);
    void add_Qtbr (size_t a, const Distribution& F, const bound_t& P);
    void add_Qee(const Distribution& F);

    // N is the Number density (inverse au^3) of particles to be added at energy e.
    void addDeltaLike(double e, double N);

    // Sets the object to have a MB distribution
    void set_maxwellian(double N, double T);

    // defines the f interpolation
    static void set_elec_points(size_t n, double min_e, double max_e);
    static size_t size;
private:
    vector<double> f;
    static BasisSet* basis;
    // static vector<double> grid; // Energy grid, Hartree
    // static double e_from_i(size_t i);
    // static size_t i_from_e(double e);
    // static double min_e, max_e; // Hartree
};

#endif /* end of include guard: RATESYSTEM_CXX_H */
