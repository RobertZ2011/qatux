#ifndef QATUX_DOUBLE_GATE
#define QATUX_DOUBLE_GATE

namespace Qatux {
    template<int Q1, int Q2, int NQUBITS, typename T>
    struct DoubleGate {
        static inline Vector<pow2(NQUBITS), T> calculate(const Vector<pow2(NQUBITS), T>& state, const Matrix<4, 4, T>& op) {
        }
    };

    template<typename T>
    struct DoubleGate<0, 1, 2, T> {
        static inline Vector<4, T> calculate(const Vector<4, T>& state, const Matrix<4, 4, T>& op) {
            return op * state;
        }
    };

    template<int NQUBITS, typename T>
    struct DoubleGate<0, 1, NQUBITS, T> {
        static inline Vector<pow2(NQUBITS), T> calculate(const Vector<pow2(NQUBITS), T>& state, const Matrix<4, 4, T>& op) {
            return kroneckerProduct(IdentityGate<NQUBITS - 2, T>::value(), op) * state;
        }
    };

    template<int NQUBITS, typename T>
    struct DoubleGate<NQUBITS - 2, NQUBITS - 1, NQUBITS, T> {
        static inline Vector<pow2(NQUBITS), T> calculate(const Vector<pow2(NQUBITS), T>& state, const Matrix<4, 4, T>& op) {
            return kroneckerProduct(op, IdentityGate<NQUBITS - 2, T>::value()) * state;
        }
    };
};
#endif
