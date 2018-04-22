#ifndef QATUX_GENERAL
#define QATUX_GENERAL

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>
#include <bitset>
#include <cmath>
#include <random>
#include <assert.h>

namespace Qatux {
    template<typename T = float>
    using Complex = std::complex<T>;

    template<typename T = float>
    using Vector = Eigen::Matrix<Complex<T>, Eigen::Dynamic, 1>;

    template<typename T = float>
    using Matrix = Eigen::SparseMatrix<Complex<T>>;

    int pow2(int exp);
    void printBin(int N, int BITS);

    //return |0>
    template<typename T = float>
    Vector<T> zero(void);

    //return |1>
    template<typename T = float>
    Vector<T> one(void);

    template<typename T = float>
    Matrix<T> identityGate(int NQUBITS);

    //returns |N> for a space of size NQUBITS
    template<typename T = float>
    Vector<T> basis(int N, int NQUBITS);
}

#endif
