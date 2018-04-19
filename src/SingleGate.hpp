#ifndef QATUX_SINGLE_GATE
#define QATUX_SINGLE_GATE

namespace Qatux {
    /*template<int Q, int NQUBITS, typename T>
    struct SingleGate {
        static inline Vector<pow2(NQUBITS), T> calculate(const Vector<pow2(NQUBITS), T>& state, const Eigen::Matrix<Complex<T>, 2, 2>& op) {
            constexpr int aboveCount = NQUBITS - (Q + 1);
            constexpr int belowCount = Q;
            return kroneckerProduct(IdentityGate<aboveCount, T>::value(), kroneckerProduct(op, IdentityGate<belowCount, T>::value())) * state;
        }
    };

    template<typename T>
    struct SingleGate<0, 1, T> {
        static inline Vector<2, T> calculate(const Vector<2, T>& state, const Eigen::Matrix<Complex<T>, 2, 2>& op) {
            return op * state;
        }
    };

    template<int NQUBITS, typename T>
    struct SingleGate<0, NQUBITS, T> {
        static inline Vector<pow2(NQUBITS), T> calculate(const Vector<pow2(NQUBITS), T>& state, const Eigen::Matrix<Complex<T>, 2, 2>& op) {
            return kroneckerProduct(IdentityGate<NQUBITS - 1, T>::value(), op) * state;
        }
    };

    template<int NQUBITS, typename T>
    struct SingleGate<NQUBITS - 1, NQUBITS, T> {
        static inline Vector<pow2(NQUBITS), T> calculate(const Vector<pow2(NQUBITS), T>& state, const Eigen::Matrix<Complex<T>, 2, 2>& op) {
            return kroneckerProduct(op, IdentityGate<NQUBITS - 1, T>::value()) * state;
        }
    };*/

    template<typename T>
    Vector<T> singleGate(int Q, int NQUBITS, const Vector<T>& state, const Matrix<T>& op) {
        if(Q == 0 && NQUBITS == 1) {
            return op * state;
        }
        else
        if(Q == NQUBITS - 1) {
            return (kroneckerProduct(op, identityGate<T>(NQUBITS - 1))) * state;
        }
        else
        if(Q == 0) {
            return (kroneckerProduct(identityGate<T>(NQUBITS - 1), op)) * state;
        }
        else {
            int aboveCount = NQUBITS - (Q + 1);
            int belowCount = Q;
            Matrix<T> higher = kroneckerProduct(identityGate<T>(aboveCount), op);
            Matrix<T> lower = identityGate<T>(belowCount);
            Matrix<T> matrix = kroneckerProduct(higher, lower);
            return  matrix * state;
        }
    }
};
#endif
