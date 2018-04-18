#ifndef QATUX_GENERAL
#define QATUX_GENERAL

#include <Eigen/KroneckerProduct>

namespace Qatux {
    template<typename T = float>
    using Complex = std::complex<T>;

    template<int N, typename T = float>
    using Vector = Eigen::Matrix<Complex<T>, N, 1>;

    template<int M, int N, typename T = float>
    using Matrix = Eigen::Matrix<Complex<T>, M, N>;

    inline constexpr int pow2(int exp) {
        return (exp == 0) ? 1 : 2 * pow2(exp - 1);
    }

    //return |0>
    template<typename T = float>
    inline Vector<2, T> zero(void) {
        Vector<2, T> result;
        result[0] = 1.0;
        result[1] = 0.0;
        return result;
    }

    //return |1>
    template<typename T = float>
    inline Vector<2, T> one(void) {
        Vector<2, T> result;
        result[0] = 0.0;
        result[1] = 1.0;
        return result;
    }

    //tensor product
    template<int N, int M, typename T>
    inline Vector<N * M, T> operator%(const Vector<N, T>& v, const Vector<M, T>& w) {
        return kroneckerProduct(v, w);
    }

    template<int NQUBITS, typename T>
    struct IdentityGate {
        static inline Matrix<pow2(NQUBITS), pow2(NQUBITS), T> value(void) {
            return kroneckerProduct(Matrix<2, 2, T>::Identity(), IdentityGate<NQUBITS - 1, T>::value());
        }
    };

    template<typename T>
    struct IdentityGate<1, T> {
        static inline Matrix<2, 2, T> value(void) {
            return Matrix<2, 2, T>::Identity();
        }
    };

    //returns |N> for a space of size NQUBITS
    template<int N, int NQUBITS, typename T = float>
    struct Basis {
        static inline Vector<pow2(NQUBITS), T> value(void) {
            static_assert(N >= 0, "Attempt to create negative qubit basis");
            static_assert(NQUBITS > 0, "Attempt to create basis for 0 or less qubit state");

            constexpr int firstBit = N >> (NQUBITS - 1) ; // N >> (NQUBITS - 1) gets the MSB of the number
            constexpr int remaining = ~(1 << (NQUBITS - 1)) & N; //gets the LSBs of the number
            return Basis<firstBit, 1, T>::value() % Basis<remaining, NQUBITS - 1, T>::value();
        }
    };


    template<typename T>
    struct Basis<0, 1, T> {
        static inline Vector<2, T> value(void) {
            return zero<T>();
        }
    };


    template<typename T>
    struct Basis<1, 1, T> {
        static inline Vector<2, T> value(void) {
            return one<T>();
        }
    };
}

#endif
