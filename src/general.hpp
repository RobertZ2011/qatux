#ifndef QATUX_GENERAL
#define QATUX_GENERAL

#include <Eigen/KroneckerProduct>

namespace Qatux {
    template<typename T = float>
    using Complex = std::complex<T>;

    template<typename T = float>
    using Vector = Eigen::Matrix<Complex<T>, Eigen::Dynamic, 1>;

    template<typename T = float>
    using Matrix = Eigen::SparseMatrix<Complex<T>>;

    inline constexpr int pow2(int exp) {
        return (exp == 0) ? 1 : 2 * pow2(exp - 1);
    }

    //return |0>
    template<typename T = float>
    inline Vector<T> zero(void) {
        Vector<T> result(2);
        result << 1.0, 0.0;
        return result;
    }

    //return |1>
    template<typename T = float>
    inline Vector<T> one(void) {
        Vector<T> result(2);
        result << 0.0, 1.0;
        return result;
    }

    template<typename T = float>
    Matrix<T> identityGate(int NQUBITS) {
        Matrix<T> gate(2, 2);
        Matrix<T> identity(2, 2);

        gate.setIdentity();
        gate.setIdentity();

        for(int i = 0; i < NQUBITS - 1; i++) {
            gate = kroneckerProduct(gate, identity).eval();
        }

        return gate;
    }

    //returns |N> for a space of size NQUBITS
    template<typename T = float>
    Vector<T> basis(int N, int NQUBITS) {
        Vector<T> basis(1);
        int highBit;
        int remaining = N;
        basis << 1.0;

        for(int i = NQUBITS - 1; i >= 0; i--) {
            highBit = remaining >> i;
            remaining = ~(1 << i) & remaining;

            if(highBit) {
                basis = kroneckerProduct(basis, one<T>()).eval();
            }
            else {
                basis = kroneckerProduct(basis, zero<T>()).eval();
            }
        }

        return basis;
    }
}

#endif
