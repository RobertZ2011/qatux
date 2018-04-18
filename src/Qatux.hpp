#ifndef QATUX
#define QATUX

#include <Eigen/Dense>
#include <bitset>
#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>

#include "general.hpp"
#include "SingleGate.hpp"

namespace Qatux {
    template<int N, typename T = float>
    class State {
    private:
        Vector<pow2(N), T> state;
        std::random_device rd;
        std::mt19937_64 generator;
        std::uniform_real_distribution<T> distribution;

    public:
        static_assert(N <= 32, "32 Qubits is the maximum supported size");

        State(const Vector<pow2(N), T>& state) {
            this->state = state;
            this->generator = std::mt19937_64(this->rd());
            this->distribution = std::uniform_real_distribution<T>(0.0, 1.0);
        }

        ~State(void) = default;

        template<int Q>
        void checkQubit(void) {
            static_assert(Q >= 0, "Negative qubit index");
            static_assert(Q <= N, "Qubit index exceeds state size");
        }

        template<int Q>
        void hadamard(void) {
            this->checkQubit<Q>();

            Eigen::Matrix<Complex<T>, 2, 2> op;
            op << 1.0,  1.0,
                  1.0, -1.0;
            op *= 1.0 / sqrt(2.0);
            this->state = SingleGate<Q, N, T>::calculate(state, op);
        }

        template<int Q>
        void notGate(void) {
            this->checkQubit<Q>();

            Eigen::Matrix<Complex<T>, 2, 2> op;
            op << 0.0, 1.0,
                  1.0, 0.0;
            this->state = SingleGate<Q, N, T>::calculate(state, op);
        }

        template<int Q>
        void phaseShift(T phi) {
            this->checkQubit<Q>();

            Eigen::Matrix<Complex<T>, 2, 2> op;
            op << 1.0, 0.0,
                  0.0, std::exp(Complex<T>(0.0, phi));

            this->state = SingleGate<Q, N, T>::calculate(state, op);
        }

        void showOutcomes(void) {
            Vector<pow2(N), T> probabilities= this->state.array() * this->state.conjugate().array();

            std::cout << "State ";

            if(N > 6) {
                std::cout << std::string(" ", N - 6);
            }

            std::cout << "| Percentage" << std::endl;

            for(int i = 0; i < pow2(N); i++) {
                std::bitset<N> state = i;
                auto flags = std::cout.flags();

                if(N < 6) {
                    std::cout << std::string(" ", 6 - N + 1);
                }

                std::cout << std::setfill('0') << std::setw(N) << state;
                std::cout.flags(flags);
                std::cout << "| " << std::fixed << std::setprecision(2) << 100 * probabilities[i].real() << std::endl;
                std::cout.flags(flags);
            }
        }

        std::bitset<N> measure(void) {
            T rand = this->distribution(this->generator);
            T start = 0;
            Vector<pow2(N), T> probabilities= this->state.array() * this->state.conjugate().array();
            int i;

            for(i = 0; i < pow2(N); i++) {
                start += probabilities[i].real();
                if(start >= rand) {
                    break;
                }
            }

            return std::bitset<N>(i);
        }
    };
}
#endif
