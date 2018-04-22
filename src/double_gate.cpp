#include "general.hpp"
#include "double_gate.hpp"

namespace Qatux {
    template<typename T>
    Matrix<T> swap2(void) {
        Matrix<T> mat(4, 4);
        mat.setIdentity();
        mat.coeffRef(1, 1) = 0.0;
        mat.coeffRef(2, 2) = 0.0;
        mat.insert(1, 2) = 1.0;
        mat.insert(2, 1) = 1.0;
        return mat;
    }

    template Matrix<float> swap2<float>(void);
    template Matrix<double> swap2<double>(void);

    template<typename T>
    Matrix<T> doubleGateConsec(int Q1, int Q2, int NQUBITS, const Matrix<T>& op) {
        if(Q1 < Q2 && NQUBITS == 2) {
            return op;
        }
        else
        if(Q1 > Q2 && NQUBITS == 2) {
            Matrix<T> swapMatrix = swap2<T>();
            return swapMatrix * op * swapMatrix;
        }
        if(Q1 == 0 && Q2 == 1) {
            return kroneckerProduct(identityGate<T>(NQUBITS - 2), op);
        }
        else
        if(Q2 == 1 && Q1 == 0) {
            return kroneckerProduct(op, identityGate<T>(NQUBITS - 2));
        }
        else {
            int lower = (Q1 < Q2) ? Q1 : Q2;
            int upperSize = NQUBITS - lower - 2;
            int lowerSize = lower;

            if(Q1 < Q2) {
                return kroneckerProduct(kroneckerProduct(identityGate<T>(upperSize), op), identityGate<T>(lowerSize));
            }
            else {
                Matrix<T> swapMatrix = swap2<T>();
                return kroneckerProduct((identityGate<T>(upperSize), swapMatrix * op * swapMatrix), identityGate<T>(lowerSize));
            }
        }
    }

    template Matrix<float> doubleGateConsec<float>(int Q1, int Q2, int NQUBITS, const Matrix<float>& op);
    template Matrix<double> doubleGateConsec<double>(int Q1, int Q2, int NQUBITS, const Matrix<double>& op);

    //generate a matrix that swaps Q1 and Q2
    template<typename T>
    Matrix<T> swapN(int Q1, int Q2, int NQUBITS) {
        Matrix<T> mat(pow2(NQUBITS), pow2(NQUBITS));
        Matrix<T> sw = swap2<T>();
        int lower = (Q1 < Q2) ? Q1 : Q2;
        int higher = (Q1 < Q2) ? Q2 : Q1;
        int above = 0;
        int below = 0;
        int requiredSwaps = higher - lower;

        mat.setIdentity();
        for(int i = 0; i < requiredSwaps; i++) {
            below = lower + i;
            above = NQUBITS - lower - i - 2;

            if(below != 0 && above != 0) {
                mat = kroneckerProduct(kroneckerProduct(identityGate<T>(above), sw), identityGate<T>(below)) * mat;
            }
            else
            if(below == 0) {
                mat = kroneckerProduct(identityGate<T>(above), sw) * mat;
            }
            else {
                mat = kroneckerProduct(sw, identityGate<T>(below)) * mat;
            }
        }

        for(int i = requiredSwaps - 2; i >= 0; i--) {
            below = lower + i;
            above = NQUBITS - lower - i - 2;

            if(below != 0 && above != 0) {
                mat = kroneckerProduct(kroneckerProduct(identityGate<T>(above), sw), identityGate<T>(below)) * mat;
            }
            else
            if(below == 0) {
                mat = kroneckerProduct(identityGate<T>(above), sw) * mat;
            }
            else {
                mat = kroneckerProduct(sw, identityGate<T>(below)) * mat;
            }
        }

        return mat;
    }

    template Matrix<float> swapN<float>(int Q1, int Q2, int NQUBITS);
    template Matrix<double> swapN<double>(int Q1, int Q2, int NQUBITS);

    template<typename T>
    Vector<T> doubleGate(int Q1, int Q2, int NQUBITS, const Vector<T>& state, const Matrix<T>& op) {
        if(Q2 == Q1 + 1 || Q1 == Q2 + 1) {
            return doubleGateConsec(Q1, Q2, NQUBITS, op) * state;
        }
        else {
            if(Q1 < Q2) {
                Matrix<T> swapMatrix = swapN<T>(Q1, Q2 - 1, NQUBITS);
                return swapMatrix * doubleGateConsec(Q2 - 1, Q2, NQUBITS, op) * swapMatrix * state;
            }
            else {
                Matrix<T> swapMatrix = swapN<T>(Q2, Q1 + 1, NQUBITS);
                auto p = doubleGateConsec(Q1, Q1 + 1, NQUBITS, op);
                return swapMatrix * p * swapMatrix * state;
            }
        }
    }

    template Vector<float> doubleGate<float>(int Q1, int Q2, int NQUBITS, const Vector<float>& state, const Matrix<float>& op);
    template Vector<double> doubleGate<double>(int Q1, int Q2, int NQUBITS, const Vector<double>& state, const Matrix<double>& op);
}
