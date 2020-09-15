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
const int GAUSS_ORDER_EII = 4;
const int GAUSS_ORDER_TBR = 4;
const double gaussX_EII[] = {-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526};
const double gaussW_EII[] = {0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538};
const double gaussX_TBR[] = {-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526};
const double gaussW_TBR[] = {0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538};


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
    SplineIntegral() {};
    void precompute_Q_coeffs(vector<RateData::Atom>& Atoms);
    // Precalculators. These delete the vectors Gamma and populate them with the calculated coefficients.)
    void Gamma_eii( eiiGraph& Gamma, const std::vector<RateData::EIIdata>& eii, size_t J) const;
    void Gamma_tbr( eiiGraph& Gamma, const std::vector<RateData::InverseEIIdata>& eii, size_t J, size_t K) const;
    // overloads that take non-vectorised inputs
    void Gamma_eii( std::vector<SparsePair>& Gamma_xi, const RateData::EIIdata& eii, size_t J) const;
    void Gamma_tbr( std::vector<SparsePair>& Gamma_xi, const RateData::InverseEIIdata& eii, size_t J, size_t K) const;

    int i_from_e(double e);
    double area(size_t idx);

    bool has_Qeii() { return _has_Qeii; }
    bool has_Qtbr() { return _has_Qtbr; }
    bool has_Qee() { return _has_Qee; }
    Q_eii_t Q_EII;
    Q_tbr_t Q_TBR;
protected:
    // Computes the overlap of the J^th df/ft term with the K^th basis function in f
    double calc_Q_eii( const RateData::EIIdata& eii, size_t J, size_t K) const;
    sparse_matrix calc_Q_tbr( const RateData::InverseEIIdata& tbr, size_t J) const;
    sparse_matrix calc_Q_ee(size_t J) const;
    bool _has_Qeii = false;  // Flags wheter Q_EII has been calculated
    bool _has_Qtbr = false;  // Flags wheter Q_TBR has been calculated
    bool _has_Qee = false;  // Flags wheter Q_EE has been calculated

    // Q_ee_t Q_EE;
    static constexpr double DBL_CUTOFF_TBR = 1e-20;
    static constexpr double DBL_CUTOFF_QEE = 1e-20;
};

#endif /* end of include guard: AC4DC_SPLINEINTEGRAL_CXX_H */
