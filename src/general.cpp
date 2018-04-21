#include "general.hpp"

namespace Qatux {
    int pow2(int exp) {
        return (exp == 0) ? 1 : 2 * pow2(exp - 1);
    }

    void printBin(int N, int BITS) {
        for(int i = BITS - 1 ; i >= 0; i--) {
            if((N >> i) & 0x1) {
                putchar('1');
            }
            else {
                putchar('0');
            }
        }
    }

    //return |0>
    template<typename T>
    Vector<T> zero(void) {
        Vector<T> result(2);
        result << 1.0, 0.0;
        return result;
    }

    template Vector<float> zero<float>(void);
    template Vector<double> zero<double>(void);

    //return |1>
    template<typename T>
    Vector<T> one(void) {
        Vector<T> result(2);
        result << 0.0, 1.0;
        return result;
    }

    template Vector<float> one<float>(void);
    template Vector<double> one<double>(void);

    template<typename T>
    Matrix<T> identityGate(int NQUBITS) {
        Matrix<T> gate(2, 2);
        Matrix<T> identity(2, 2);

        gate.setIdentity();
        identity.setIdentity();

        for(int i = 0; i < NQUBITS - 1; i++) {
            gate = kroneckerProduct(gate, identity).eval();
        }

        return gate;
    }

    template Matrix<float> identityGate<float>(int NQUBITS);
    template Matrix<double> identityGate<double>(int NQUBITS);

    //returns |N> for a space of size NQUBITS
    template<typename T>
    Vector<T> basis(int N, int NQUBITS) {
        Vector<T> basis(1);
        int highBit;
        int remaining = N;
        basis << 1.0;

        for(int i = NQUBITS - 1; i >= 0; i--) {
            highBit = remaining >> i;
            remaining = ~(1 << i) & remaining;

            if(highBit) {
                basis = kroneckerProduct(basis, one<T>()).eval();
            }
            else {
                basis = kroneckerProduct(basis, zero<T>()).eval();
            }
        }

        return basis;
    }

    template Vector<float> basis<float>(int N, int NQUBITS);
    template Vector<double> basis<double>(int N, int NQUBITS);
}
