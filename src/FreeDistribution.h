#ifndef FREEDISTRIBUTION_CXX_H
#define FREEDISTRIBUTION_CXX_H


#include <sstream>
#include <assert.h>
#include <iostream>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SparseCholesky>
#include "Constant.h"

#define BSPLINE_ORDER 4 // Cubics should be fine, increase if needed
#define GAUSS_ORDER 10

class BasisSet{
public:
    BasisSet(){};
    void set_parameters(size_t n, double min, double max);
    Eigen::VectorXd Sinv(Eigen::VectorXd deltaf);
    double operator()(size_t i, double x);
    inline double supp_max(unsigned k);
    inline double supp_min(unsigned k);
    size_t gridlen(){
        return knot.size();
    }
    double grid(size_t i){
        return knot[i];
    }
    size_t num_funcs;
private:
    // Eigen::LDLT<Eigen::MatrixXd > cholmachine;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > cholmachine;
    std::vector<double> knot;
    double overlap(size_t j, size_t k);
};

// Represents a statistical distribution of electrons. Internal units are atomic units.
class Distribution
{
public:
    Distribution() {
        f.resize(size);
    }

    inline double operator[](size_t n) const{
        return this->f[n];
    }

    // vector-space algebra
    Distribution& operator+=(const Distribution& d){
        for (size_t i=0; i<size; i++) {
            f[i] += d.f[i];
        }
        return *this;
    }

    Distribution& operator*=(double x){
        for (int i=0; i<f.size(); i++){
            f[i] *= x;
        }
        return *this;
    }

    Distribution& operator=(const Distribution& d){
        f = d.f;
        return *this;
    }

    Distribution& operator=(double y){
        for (int i=0; i<f.size(); i++){
            f[i] = y;
        }
        return *this;
    }

    static void Gamma_eii(Eigen::SparseMatrix<double>& Gamma, const CustomDataType::EIIdata& eii, size_t K, size_t a);
    static void Gamma_tbr(Eigen::SparseMatrix<double>& Gamma, const CustomDataType::EIIdata& eii, size_t J, size_t K, size_t a);
    void add_Qeii (size_t a, const Distribution& F, const bound_t& P);
    void add_Qtbr (size_t a, const Distribution& F, const bound_t& P);
    void add_Qee(const Distribution& F);

    // N is the Number density (inverse au^3) of particles to be added at energy e.
    void addDeltaLike(double e, double N);

    // Sets the object to have a MB distribution
    void set_maxwellian(double N, double T);

    static std::string get_energies_eV();
    // defines the f interpolation
    static void set_elec_points(size_t n, double min_e, double max_e);
    static size_t size;
private:
    // Eigen::VectorXd f;
    std::vector<double> f;
    static BasisSet basis;
    // static vector<double> grid; // Energy grid, Hartree
    // static double e_from_i(size_t i);
    // static size_t i_from_e(double e);
    // static double min_e, max_e; // Hartree
};

#endif /* end of include guard: RATESYSTEM_CXX_H */
