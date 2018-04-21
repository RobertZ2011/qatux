#ifndef QATUX_SINGLE_GATE
#define QATUX_SINGLE_GATE

namespace Qatux {
    template<typename T>
    Vector<T> singleGate(int Q, int NQUBITS, const Vector<T>& state, const Matrix<T>& op);
};
#endif
