#ifndef QATUX_SINGLE_GATE
#define QATUX_SINGLE_GATE

namespace Qatux {
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
