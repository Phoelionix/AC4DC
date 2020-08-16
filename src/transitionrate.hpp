#ifndef AS_BASICMATRIX_CXX_H
#define AS_BASICMATRIX_CXX_H
#include <vector>
#include <assert.h>

typedef std::vector<double> bound_t; // Probabilities of state

namespace {
namespace GammaType{
    class TransitionRate{
    public:
        TransitionRate(){}; // makes an empty matrix
        TransitionRate(size_t n);
        void resize(size_t n);
        double& operator()(size_t i, size_t j);
        // double at(size_ti, size_t j);
        bound_t operator*(const bound_t& v) const;
        bound_t calc_delta(const bound_t& P) const;
        TransitionRate& operator=(const TransitionRate& M);
        void print();
    private:
        std::vector<bound_t> _matrix;
        size_t N;

    };


    void TransitionRate::resize(size_t n){
        N = n;

        _matrix.resize(n);
        for (auto& v : _matrix){
            v.resize(n);
        }
    }


    TransitionRate::TransitionRate(size_t n){
        resize(n);
    }

    double& TransitionRate::operator()(size_t i, size_t j){
        return _matrix[i][j];
    }

    bound_t TransitionRate::operator*(const bound_t& v) const {
        // A_ij v_j
        bound_t retval(N);
        for (int i=0; i<N; i++){
            retval[i] = 0;
            for (int j=0; j<N; j++){
                retval[i] += _matrix[i][j]*v[j];
            }
        }
        return retval;
    }

    TransitionRate& TransitionRate::operator=(const TransitionRate& M) {
        for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                _matrix[i][j] = M._matrix[i][j];
            }
        }
        return *this;
    }

    bound_t TransitionRate::calc_delta(const bound_t& P) const{
        assert(P.size() == N);
        bound_t Pdot(N);
        // Compute the changes in P
        for (size_t i = 0; i < N; i++) {

            for (size_t j = 0; j < i; j++) {
                Pdot[i] += _matrix[i][j] * P[j] - _matrix[j][i] * P[i];
            }
            // Avoid j=i
            for (size_t j = i+1; j < N; j++) {
                Pdot[i] += _matrix[i][j] * P[j] - _matrix[j][i] * P[i];
            }
        }
        return Pdot;
    }
};
};

#endif /* end of include guard: AS_BASICMATRIX_CXX_H */
