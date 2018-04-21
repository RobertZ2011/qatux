#ifndef QATUX
#define QATUX

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <bitset>
#include <cmath>
#include <random>
#include <assert.h>

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
        State(int NQUBITS, const Vector<T>& state) : state(state), NQUBITS(NQUBITS){
            assert(state.rows() == pow2(NQUBITS));
            this->generator = std::mt19937_64(this->rd());
            this->distribution = std::uniform_real_distribution<T>(0.0, 1.0);
        }

        ~State(void) = default;

        const State& operator=(const Vector<T>& state) {
            this->state = state;
            return *this;
        }

        void hadamard(int Q) {
            Matrix<T> op(2, 2);

            op.insert(0, 0) = 1.0;
            op.insert(0, 1) = 1.0;
            op.insert(1, 0) = 1.0;
            op.insert(1, 1) = -1.0;

            op *= 1.0 / sqrt(2.0);
            this->state = singleGate(Q, this->NQUBITS, this->state, op);
        }

        void notGate(int Q) {
            Matrix<T> op(2, 2);

            op.insert(0, 0) = 0.0;
            op.insert(0, 1) = 1.0;
            op.insert(1, 0) = 1.0;
            op.insert(1, 1) = 0.0;

            this->state = singleGate(Q, this->NQUBITS, this->state, op);
        }

        void phaseShift(T phi) {
            Matrix<T> op(2, 2);

            op.insert(0, 0) = 1.0;
            op.insert(0, 1) = 0.0;
            op.insert(1, 0) = 0.0;
            op.insert(1, 1) = std::exp(Complex<T>(0.0, phi));

            this->state = singleGate(state, op);
        }

        void controlled(int Q1, int Q2, const Matrix<T>& base) {
            Matrix<T> op(4, 4);

            op.reserve(6);
            op.setIdentity();

            op.insert(2, 3) = base.coeff(0, 1);
            op.insert(3, 2) = base.coeff(1, 0);

            op.coeffRef(2, 2) = base.coeff(0, 0);
            op.coeffRef(3, 3) = base.coeff(1, 1);

            this->state = doubleGate<T>(Q1, Q2, this->NQUBITS, state, op);
        }

        void cnot(int Q1, int Q2) {
            Matrix<T> op(2, 2);

            op.reserve(2);

            op.insert(0, 0) = 0.0;
            op.insert(0, 1) = 1.0;
            op.insert(1, 0) = 1.0;
            op.insert(1, 1) = 0.0;

            this->controlled(Q1, Q2, op);
        }

        void swap(int Q1, int Q2) {
            if(this->NQUBITS == 2) {
                this->state = swap2<T>() * this->state;
            }
            else
            if(Q2 == Q1 + 1 || Q1 == Q2 + 1) {
                this->state = doubleGateConsec(Q1, Q2, NQUBITS, swap2<T>());
            }
            else {
                this->state = swapN<T>(Q1, Q2, this->NQUBITS) * this->state;
            }
        }

        void showOutcomes(void) {
            Vector<T> probabilities  = this->state.array() * this->state.conjugate().array();

            printf("State | Percentage\n");
            for(int i = 0; i < this->state.rows(); i++) {
                printBin(i, this->NQUBITS);
                printf(" | %.2f\n", (float) probabilities[i].real() * 100.0);
            }
        }

        std::bitset<32> measure(void) {
            T rand = this->distribution(this->generator);
            T start = 0;
            Vector<T> probabilities= this->state.array() * this->state.conjugate().array();
            int i;

            for(i = 0; i < pow2(this->NQUBITS); i++) {
                start += probabilities[i].real();
                if(start >= rand) {
                    break;
                }
            }

            return std::bitset<32>(i);
        }
    };
}
#endif
