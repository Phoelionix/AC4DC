/**
 * @file SplineIntegral.h
 * @brief Part of Sanders' continuum plasma extension.
 * 
 */

/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================*/

#ifndef AC4DC_SPLINEINTEGRAL_CXX_H
#define AC4DC_SPLINEINTEGRAL_CXX_H

#include "Constant.h"
#include "Dipole.h"
#include "SplineBasis.h"
#include <omp.h>

struct SparsePair
{
    int idx;
    double val;
    SparsePair() {idx=0; val=0;};
    SparsePair( int i, double v ) : idx(i), val(v) {};
    ~SparsePair() {};
    SparsePair& operator=(SparsePair tup) {
        idx = tup.idx;
        val = tup.val;
        return *this;
    }
};

// defines for numerical integration
static const int GAUSS_ORDER_EII = 10;
// static const double gaussX_EII[] = {-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526};
// static const double gaussW_EII[] = {0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538};
static const int GAUSS_ORDER_TBR = 10;
// static const double gaussX_TBR[] = {-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526};
// static const double gaussW_TBR[] = {0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538};
static const double gaussX_EE[] = {-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526};
static const double gaussW_EE[] = {0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538};
static const int GAUSS_ORDER_EE = 4;
// #define gaussX_EE gaussX_10
// #define gaussW_EE gaussW_10
#define gaussX_EII gaussX_10
#define gaussW_EII gaussW_10
#define gaussX_TBR gaussX_10
#define gaussW_TBR gaussW_10

// Specialised quadrature rules for integrals of the form \int_-1^1 dx sqrt(x+1) f(x)
// calculated using this excelent tool https://keisan.casio.com/exec/system/15041456505628
static const int GAUSS_ORDER_SQRT = 6; // i.e. length of the below containers -S.P.
static const double gaussX_sqrt[] = {-0.8937779292142465272226, -0.5977085045355194007473, -0.1747746522411171589491, 0.2850548710873168381781, 0.6839736445131432976427, 0.9372325703904229510983};
static const double gaussW_sqrt[] = {0.0679848323648580100291, 0.2364639422771327772543, 0.4158087090592895378275, 0.5047613331425565232555, 0.4387744028120991078995, 0.2218248635081907754697};

// Same as above but for things that resemble 1/sqrt(e)
static const int GAUSS_ORDER_ISQRT = 6;
static const double gaussX_isqrt[] = {-0.9686331867851990899941, -0.7293999766895035262393, -0.3101152411451654769277, 0.1855002554630833445573, 0.634856026533749915623, 0.926922557405643528633};
static const double gaussW_isqrt[] = {0.704694262429010138266, 0.6604166237708736260642, 0.5746442606084036114914, 0.4527698865360327398879, 0.302470090347234013247, 0.1334320010546359686468};

// Interpretation:
// eiiGraph[xi] -> vector of transitions away from configuration xi, e.g.
// eiiGraph[1] = [ (3, 0.33), (4,0.22)] 
// i.e. there are two transitions away from config 1, towards configs 3 and 4
// each corresponding to gamma values of 0.33 and 0.22
// (these must still be multiplied by the free distribution)
typedef std::vector<std::vector<SparsePair> > eiiGraph;

struct SparseTriple
{
    int K;
    int L;
    double val;
    static SparseTriple Zero() {
        SparseTriple tmp;
        tmp.K=0;
        tmp.L=0;
        tmp.val=0;
        return tmp;
    }
};

class SplineIntegral : public BasisSet {
public:
    typedef std::vector<SparseTriple> sparse_matrix;
    typedef std::vector<SparsePair> pair_list;

    typedef std::vector<std::vector< std::vector< std::vector<double> > > > Q_eii_t;
    // Interpretation: J^th matrix element of dQ/dt given by
    // Q_eii[a][xi][J][K] * P^a[xi] * F[K]
    //       ^  ^   ^  ^
    //       |  |   Basis indices (J index for LHS, K index for RHS)
    //       |  state within atom a
    //       Atom
    typedef std::vector<std::vector<std::vector<sparse_matrix> > > Q_tbr_t;
    // Interpretation: J^th matrix element of dQ/dt given by
    // vector<SparseTriple> &nv = Q_tbr[a][xi][J]
    // for (auto& Q : nv) {
    //      tmp += nv.val * P^a[xi] * F[nv.K] * F[nv.L]
    // }
    typedef std::vector<std::vector<pair_list> > Q_ee_t;
    // Interpretation: J'th element of df/dt is given by
    // for (auto& q : Q_EE[J][K]) {
    //     tmp += q.val*F[K]*F[q.idx]
    // }
    // 1000 times fewer components than QTBR, not a problem

     /**
     * @brief // Replaces inner knots with the ones provided, using the current boundaries to determine which are inner. 
     * @details Note that knot[i] is just an energy. Intended for use in loading prev. simulation states.
     * 
     * @param inner_knots Inner knot energies in Hartrees [Ha]. The new (non-boundary) knot vector, in case you knew not.
     * @return SplineIntegral& 
     */
    SplineIntegral& operator=(vector<double> inner_knots){
        // Count number of boundary knots
        int lower = 0;
        int upper = 0;
        for(size_t k = 0; k  < knot.size(); k++){
            if(knot[k] <= _min){  
                lower++;}
            else if (knot[k] > _max){ // note knot at _max is not boundary knot. knot.
                upper++;}
        }
        // Remove inner knots
        knot.erase(std::next(knot.begin(), lower), std::next(knot.begin(), knot.size()- upper));
        // Insert new knots
        knot.reserve(lower + inner_knots.size() + upper);
        knot.insert(knot.begin() + lower,inner_knots.begin(),inner_knots.end());
        return *this;        
    }    


    vector<double> get_knot(){return knot;}

    SplineIntegral() {};
    void precompute_QEII_coeffs(vector<RateData::Atom>& Atoms);
    void precompute_QTBR_coeffs(vector<RateData::Atom>& Atoms);
    void precompute_QEE_coeffs();
    // Precalculators. These delete the vectors Gamma and populate them with the calculated coefficients.)
    void Gamma_eii( eiiGraph& Gamma, const std::vector<RateData::EIIdata>& eii, size_t J) const;
    void Gamma_tbr( eiiGraph& Gamma, const std::vector<RateData::InverseEIIdata>& eii, size_t J, size_t K) const;
    // overloads that take non-vectorised inputs
    void Gamma_eii( std::vector<SparsePair>& Gamma_xi, const RateData::EIIdata& eii, size_t J) const;
    void Gamma_tbr( std::vector<SparsePair>& Gamma_xi, const RateData::InverseEIIdata& eii, size_t J, size_t K) const;

    bool has_Qeii() { return _has_Qeii; }
    bool has_Qtbr() { return _has_Qtbr; }
    bool has_Qee() { return _has_Qee; }
    Q_eii_t Q_EII;
    Q_tbr_t Q_TBR;
    Q_ee_t Q_EE;
protected:
    // Computes the overlap of the J^th df/ft term with the K^th basis function in f
    double calc_Q_eii( const RateData::EIIdata& eii, size_t J, size_t K) const;
    sparse_matrix calc_Q_tbr( const RateData::InverseEIIdata& tbr, size_t J) const;
    pair_list calc_Q_ee(size_t J, size_t K) const;

    inline double Q_ee_F(double e, size_t K) const;
    inline double Q_ee_G(double e, size_t K) const;
    // sparse_matrix calc_Q_ee(size_t J) const;
    bool _has_Qeii = false;  // Flags wheter Q_EII has been calculated
    bool _has_Qtbr = false;  // Flags wheter Q_TBR has been calculated
    bool _has_Qee = false;  // Flags wheter Q_EE has been calculated

    static constexpr double DBL_CUTOFF_TBR = 1e-16; // 1e-16 is smallest machine pracision -S.P.
    static constexpr double DBL_CUTOFF_QEE = 1e-16;
};

#endif /* end of include guard: AC4DC_SPLINEINTEGRAL_CXX_H */
