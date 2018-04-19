#ifndef QATUX
#define QATUX

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <bitset>
#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>

#include "general.hpp"
#include "SingleGate.hpp"

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
        State(int NQUBITS) {
            this->state = Vector<T>(pow2(NQUBITS));
            this->NQUBITS = NQUBITS;
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

            /*op << 1.0,  1.0,
                  1.0, -1.0;*/
            op *= 1.0 / sqrt(2.0);
            this->state = singleGate(Q, this->NQUBITS, this->state, op);
        }

        /*void notGate(int Q) {
            Matrix<T> op(2, 2);
            op << 0.0, 1.0,
                  1.0, 0.0;
            this->state = singleGate(Q, this->NQUBITS, this->state, op);
        }*

        /*template<int Q>
        void phaseShift(T phi) {
            this->checkQubit<Q>();

            Eigen::Matrix<Complex<T>, 2, 2> op;
            op << 1.0, 0.0,
                  0.0, std::exp(Complex<T>(0.0, phi));

            this->state = SingleGate<Q, N, T>::calculate(state, op);
        }

        template<int Q1, int Q2>
        void controlled(const Matrix<2, 2, T>& base) {
            Matrix<4, 4, T> op = Matrix<4, 4, T>::Identity();
            op.bottomRightCorner<2, 2>() = base;

            this->state = DoubleGate<Q1, Q2, N, T>::calculate(state, op);
        }

        template<int Q1, int Q2>
        void cnot(void) {
            this->checkQubit2<Q1, Q2>();

            Matrix<2, 2, T> op;
            op << 0.0, 1.0,
                  1.0, 0.0;

            this->controlled<Q1, Q2>(op);
        }*/

        void showOutcomes(void) {
            Vector<T> probabilities= this->state.array() * this->state.conjugate().array();

            std::cout << "State ";

            /*if(N > 6) {
                std::cout << std::string(" ", N - 6);
            }*/

            std::cout << "| Percentage" << std::endl;

            for(int i = 0; i < this->state.rows(); i++) {
                std::bitset<3> state = i;
                auto flags = std::cout.flags();

                /*if(N < 6) {
                    std::cout << std::string(" ", 6 - N + 1);
                }*/

                std::cout << std::setfill('0') << std::setw(8) << state;
                std::cout.flags(flags);
                std::cout << "| " << std::fixed << std::setprecision(2) << 100 * probabilities[i].real() << std::endl;
                std::cout.flags(flags);
            }
        }

        /*std::bitset<N> measure(void) {
            T rand = this->distribution(this->generator);
            T start = 0;
            Vector<T> probabilities= this->state.array() * this->state.conjugate().array();
            int i;

            for(i = 0; i < pow2(N); i++) {
                start += probabilities[i].real();
                if(start >= rand) {
                    break;
                }
            }

            return std::bitset<N>(i);
        }*/
    };
}
#endif
