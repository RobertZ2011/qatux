#ifndef QATUX_DOUBLE_GATE
#define QATUX_DOUBLE_GATE

namespace Qatux {
    template<typename T = float>
    Matrix<T> swap2(void);

    template<typename T>
    Matrix<T> doubleGateConsec(int Q1, int Q2, int NQUBITS, const Matrix<T>& op);

    //generate a matrix that swaps Q1 and Q2
    template<typename T>
    Matrix<T> swapN(int Q1, int Q2, int NQUBITS);

    template<typename T>
    Vector<T> doubleGate(int Q1, int Q2, int NQUBITS, const Vector<T>& state, const Matrix<T>& op);
};
#endif
