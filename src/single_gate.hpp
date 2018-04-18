#ifndef QATUX_SINGLE_GATE
#define QATUX_SINGLE_GATE

namespace Qatux {
    template<int Q, int NQUBITS, typename T>
    struct SingleGateOp {
        static inline Vector<pow2(NQUBITS), T> calculate(const Vector<pow2(NQUBITS), T>& state, const Eigen::Matrix<Complex<T>, 2, 2>& op) {
            constexpr int aboveCount = NQUBITS - (Q + 1);
            constexpr int belowCount = Q;
            return kroneckerProduct(IdentityGate<aboveCount, T>::value(), kroneckerProduct(op, IdentityGate<belowCount, T>::value())) * state;
        }
    };

    template<typename T>
    struct SingleGateOp<0, 1, T> {
        static inline Vector<2, T> calculate(const Vector<2, T>& state, const Eigen::Matrix<Complex<T>, 2, 2>& op) {
            return op * state;
        }
    };

    template<int NQUBITS, typename T>
    struct SingleGateOp<0, NQUBITS, T> {
        static inline Vector<pow2(NQUBITS), T> calculate(const Vector<pow2(NQUBITS), T>& state, const Eigen::Matrix<Complex<T>, 2, 2>& op) {
            return kroneckerProduct(IdentityGate<NQUBITS - 1, T>::value(), op) * state;
        }
    };

    template<typename T>
    struct SingleGateOp<1, 2, T> {
        static inline Vector<4, T> calculate(const Vector<4, T>& state, const Eigen::Matrix<Complex<T>, 2, 2>& op) {
            return kroneckerProduct(op, IdentityGate<1, T>::value()) * state;
        }
    };

    template<int NQUBITS, typename T>
    struct SingleGateOp<NQUBITS - 1, NQUBITS, T> {
        static inline Vector<pow2(NQUBITS), T> calculate(const Vector<pow2(NQUBITS), T>& state, const Eigen::Matrix<Complex<T>, 2, 2>& op) {
            return kroneckerProduct(op, IdentityGate<NQUBITS - 1, T>::value()) * state;
        }
    };
};
#endif
