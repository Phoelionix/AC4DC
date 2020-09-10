#ifndef FREEDISTRIBUTION_CXX_H
#define FREEDISTRIBUTION_CXX_H


#include <sstream>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SparseCholesky>
#include "Constant.h"
#include "SplineIntegral.h"
#include "Dipole.h"


// Represents a statistical distribution of electrons. Internal units are atomic units.
class Distribution
{
public:
    Distribution() {
        f.resize(size);
    }

    inline double& operator[](size_t n){
        return this->f[n];
    }

    inline double operator[](size_t n) const{
        return this->f[n];
    }

    // vector-space algebra
    Distribution& operator+=(const Distribution& d){
        for (size_t i=0; i<size; i++) {
            f[i] += d.f[i];
        }
        // total += d.total;
        return *this;
    }

    Distribution& operator*=(double x){
        for (int i=0; i<f.size(); i++){
            f[i] *= x;
        }
        // total *= d.total;
        return *this;
    }

    Distribution& operator=(const Distribution& d){
        f = d.f;
        // total = d.total;
        return *this;
    }

    Distribution& operator=(double y){
        // total = y;
        for (int i=0; i<f.size(); i++){
            f[i] = y;
        }
        return *this;
    }

    double norm() const;

    // modifiers

    // applies the df/dt vector v to the overall distribution
    void applyDelta(const Eigen::VectorXd& dfdt);
    

    // Q functions
    // Computes the dfdt vector v based on internal f
    // e.g. dfdt v; F.calc_Qee(v);
    void get_Q_eii (Eigen::VectorXd& v, size_t a, const bound_t& P) const;
    void get_Q_tbr (Eigen::VectorXd& v, size_t a, const bound_t& P) const;
    void apply_Qee  (Eigen::VectorXd& v) const;
    // N is the Number density (inverse au^3) of particles to be added at energy e.
    static void addDeltaLike(Eigen::VectorXd& v, double e, double N);
    // Adds a Dirac delta to the distribution
    void addDeltaSpike(double N, double e);
    // Sets the object to have a MB distribution
    void set_maxwellian(double N, double T);

    // Precalculators
    static void Gamma_eii( eiiGraph& Gamma, const std::vector<RateData::EIIdata>& eii, size_t J){
        return basis.Gamma_eii(Gamma, eii, J);
    }
    static void Gamma_tbr( eiiGraph& Gamma, const std::vector<RateData::EIIdata>& eii, size_t J, size_t K){
        return basis.Gamma_tbr(Gamma, eii, J, K);
    }
    static void precompute_Q_coeffs(vector<RateData::Atom>& Store){
        basis.precompute_Q_coeffs(Store);
    }

    double integral(double (f)(double));

    static std::string output_energies_eV(size_t num_pts);
    std::string output_densities(size_t num_pts) const;

    double operator()(double e) const;

    // The setup function
    static void set_elec_points(size_t n, double min_e, double max_e, GridSpacing grid_style);



    static size_t size;
    static SplineIntegral basis;
private:
    // double total;
    std::vector<double> f;
    // static SplineIntegral basis;

};

ostream& operator<<(ostream& os, const Distribution& dist);

#endif /* end of include guard: RATESYSTEM_CXX_H */
