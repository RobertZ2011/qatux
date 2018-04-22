#include "general.hpp"
#include "single_gate.hpp"

namespace Qatux {
    template<typename T>
    Vector<T> singleGate(int Q, int NQUBITS, const Vector<T>& state, const Matrix<T>& op) {
        Matrix<T> mat;

        if(Q == 0 && NQUBITS == 1) {
            mat = op;
        }
        else
        if(Q == NQUBITS - 1) {
            mat = kroneckerProduct(op, identityGate<T>(NQUBITS - 1));
        }
        else
        if(Q == 0) {
            mat = kroneckerProduct(identityGate<T>(NQUBITS - 1), op);
        }
        else {
            int aboveCount = NQUBITS - (Q + 1);
            int belowCount = Q;
            Matrix<T> higher = kroneckerProduct(identityGate<T>(aboveCount), op);
            Matrix<T> lower = identityGate<T>(belowCount);
            mat = kroneckerProduct(higher, lower);
        }

        return mat * state;
    }

    template Vector<float> singleGate<float>(int Q, int NQUBITS, const Vector<float>& state, const Matrix<float>& op);
    template Vector<double> singleGate<double>(int Q, int NQUBITS, const Vector<double>& state, const Matrix<double>& op);
}
