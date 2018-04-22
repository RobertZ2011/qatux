#ifndef QATUX
#define QATUX

#include "general.hpp"
#include "single_gate.hpp"
#include "double_gate.hpp"

namespace Qatux {
    template<typename T = float>
    class State {
    private:
        Vector<T> state;
        int NQUBITS;
        std::random_device rd;
        std::mt19937_64 generator;
        std::uniform_real_distribution<T> distribution;

    public:
        State(int NQUBITS, const Vector<T>& state);
        ~State(void) = default;

        const State& operator=(const Vector<T>& state);
        void hadamard(int Q);
        void notGate(int Q);
        void phaseShift(int Q, T phi);
        void controlled(int Q1, int Q2, const Matrix<T>& base);
        void cnot(int Q1, int Q2);
        void swap(int Q1, int Q2);

        void showOutcomes(void);
        std::bitset<32> measure(void);
    };
}
#endif
