#ifndef FREEDISTRIBUTION_CXX_H
#define FREEDISTRIBUTION_CXX_H


#include <sstream>
#include <assert.h>
#include <iostream>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SparseCholesky>
#include "Constant.h"
#include "SplineBasis.h"
#include "Constant.h"
#include "Dipole.h"
#include "FreeDistribution.h"
#include "SplineBasis.h"
#include "transitionrate.hpp"

#define DBL_CUTOFF_TBR 1e-20

struct SparseEntry
{
    int K;
    int L;
    double val;
    static SparseEntry Zero(){
        SparseEntry tmp;
        tmp.K=0;
        tmp.L=0;
        tmp.val=0;
        return tmp;
    }
};

// Represents a statistical distribution of electrons. Internal units are atomic units.
class Distribution
{
public:
    typedef std::vector<SparseEntry> tbr_sparse;

    typedef std::vector<std::vector< std::vector< std::vector<double> > > > Q_eii_t;
    // Interpretation: J^th matrix element of dQ/dt given by
    // Q_eii[a][xi][J][K] * P^a[xi] * F[K]
    //       ^  ^   ^  ^
    //       |  |   Basis indices (J index for LHS, K index for RHS)
    //       |  state within atom a
    //       Atom
    typedef std::vector<std::vector<std::vector<tbr_sparse> > > Q_tbr_t;
    // Interpretation: J^th matrix element of dQ/dt given by
    // vector<SparseEntry> &nv = Q_tbr[a][xi][J]
    // for (auto& Q : nv){
    //      tmp += nv.val * P^a[xi] * F[nv.K] * F[nv.L]
    // }

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

    // Precalculators
    static void Gamma_eii( GammaType::eiiGraph& Gamma, const std::vector<RateData::EIIdata>& eii, size_t J);
    static void Gamma_tbr( GammaType::eiiGraph& Gamma, const std::vector<RateData::EIIdata>& eii, size_t J, size_t K);

    // changes v, doues NOT
    // N is the Number density (inverse au^3) of particles to be added at energy e.
    static void addDeltaLike(Eigen::VectorXd& v, double e, double N);

    static std::string get_energies_eV();

    // The setup function
    static void set_elec_points(size_t n, double min_e, double max_e);
    static void precompute_Q_coeffs(vector<RateData::Atom>& Store);
    static size_t size;
private:
    // double total;
    std::vector<double> f;
    static BasisSet basis;

    static bool has_Qeii;  // Flags wheter Q_EII has been calculated
    static bool has_Qtbr;  // Flags wheter Q_TBR has been calculated
    static Q_eii_t Q_EII;
    static Q_tbr_t Q_TBR;

    static double calc_Q_eii( const RateData::EIIdata& eii, size_t J, size_t K);
    static tbr_sparse calc_Q_tbr( const RateData::InverseEIIdata& tbr, size_t J);

};

#endif /* end of include guard: RATESYSTEM_CXX_H */
