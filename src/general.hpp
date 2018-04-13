#ifndef QATUX_GENERAL
#define QATUX_GENERAL

namespace Qatux {
    template<int N, typename T>
    using Vector = Eigen::Matrix<std::complex<T>, N, 1>;

    inline constexpr int pow2(int exp) {
        return (exp == 0) ? 1 : 2 * pow2(exp - 1);
    }

    template<typename T>
    Vector<2, T> zero(void) {
        Vector<2, T> result;
        result[0] = std::complex<T>(1.0, 0.0);
        result[1] = std::complex<T>(0.0, 0.0);
        return result;
    }

    template<typename T>
    Vector<2, T> one(void) {
        Vector<2, T> result;
        result[0] = std::complex<T>(0.0, 0.0);
        result[1] = std::complex<T>(1.0, 0.0);
        return result;
    }

    template<typename T>
    Vector<2, T> supPlus(void) {
        T norm = 1.0 / sqrt(2.0);
        Vector<2, T> result = norm * (zero<T>() + one<T>());
        return result;
    }

    template<typename T>
    Vector<2, T> supMinus(void) {
        T norm = 1.0 / sqrt(2.0);
        Vector<2, T> result = norm * (zero<T>() + one<T>());
        return result;
    }


    //tensor product
    template<int N, int M, typename T>
    Vector<N * M, T> operator%(const Vector<N, T>& n, const Vector<M, T>& m) {
        Vector<N * M, T> result;

        for(int i = 0; i < N; i++) {
            for(int j = 0; j < M; j++) {
                result[2 * i + j] = n[i] * m[j];
            }
        }

        return result;
    }

    template<int N, int NQUBITS, typename T>
    struct basis {
        static Vector<pow2(NQUBITS), T> value(void) {
            static_assert(N >= 0, "Attempt to create negative qubit basis");
            static_assert(NQUBITS > 0, "Attempt to create basis for 0 or less qubit state");

            constexpr int firstBit = N >> (NQUBITS - 1) ; // N >> (NQUBITS - 1) gets the MSB of the number
            constexpr int remaining = ~(1 << (NQUBITS - 1)) & N; //gets the LSBs of the number
            return basis<remaining, NQUBITS - 1, T>::value() % basis<firstBit, 1, T>::value();
        }
    };


    template<typename T>
    struct basis<0, 1, T> {
        static Vector<2, T> value(void) {
            return zero<T>();
        }
    };


    template<typename T>
    struct basis<1, 1, T> {
        static Vector<2, T> value(void) {
            return one<T>();
        }
    };
}

#endif
