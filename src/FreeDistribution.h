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
    // Sets the object to have a MB distribution
    void set_maxwellian(double N, double T);

    // Q functions
    // Computes the dfdt vector v based on internal f
    // e.g. dfdt v; F.calc_Qee(v);
    void get_Q_eii (Eigen::VectorXd& v, size_t a, const bound_t& P) const;
    void get_Q_tbr (Eigen::VectorXd& v, size_t a, const bound_t& P) const;
    void apply_Qee  (Eigen::VectorXd& v) const;
    // N is the Number density (inverse au^3) of particles to be added at energy e.
    static void addDeltaLike(Eigen::VectorXd& v, double e, double N);
    void addDeltaSpike(double e, double N);

    // Precalculators
    static void Gamma_eii( eiiGraph& Gamma, const std::vector<RateData::EIIdata>& eii, size_t J);
    static void Gamma_tbr( eiiGraph& Gamma, const std::vector<RateData::EIIdata>& eii, size_t J, size_t K);

    static std::string get_energies_eV();

    // The setup function
    static void set_elec_points(size_t n, double min_e, double max_e);
    static void precompute_Q_coeffs(vector<RateData::Atom>& Store);
    static size_t size;
private:
    // double total;
    std::vector<double> f;
    static SplineIntegral basis;

};

#endif /* end of include guard: RATESYSTEM_CXX_H */
