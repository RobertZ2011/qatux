#ifndef QATUX_SINGLE_GATE
#define QATUX_SINGLE_GATE

namespace Qatux {
    template<int NQUBITS, typename T>
    using SingleGateCoeffVector = Eigen::Array<Vector<2, T>, pow2(NQUBITS - 1), 1>;

    template<int Q, int NQUBITS, typename T>
    struct SingleGateCoeffVec {
            static SingleGateCoeffVector<NQUBITS, T> calculate(const Vector<pow2(NQUBITS), T>& state, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
                SingleGateCoeffVector<NQUBITS, T> coeffs;

                for(int i = 0; i < pow2(NQUBITS); i += 2) {
                    coeffs[i / 2] = state[i] * zero<T>() + state[i + 1] * one<T>();
                }

                return coeffs;
            }
    };

    template<int NQUBITS, typename T>
    struct SingleGateCoeffVec<0, NQUBITS, T> {
        static SingleGateCoeffVector<NQUBITS, T> calculate(const Vector<pow2(NQUBITS), T>& state, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            SingleGateCoeffVector<NQUBITS, T> coeffs;

            for(int i = 0; i < pow2(NQUBITS); i += 2) {
                coeffs[i / 2] = state[i] * Gzero + state[i + 1] * Gone;
            }

            return coeffs;
        }
    };

    template<int NQUBITS, typename T>
    struct SingleGateCoeffVec<NQUBITS - 1, NQUBITS, T> {
        static SingleGateCoeffVector<NQUBITS, T> calculate(const Vector<pow2(NQUBITS), T>& state, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            SingleGateCoeffVector<NQUBITS, T> coeffs;
            Vector<pow2(NQUBITS - 1), T> low = state.template segment<pow2(NQUBITS - 1)>(0);
            Vector<pow2(NQUBITS - 1), T> high = state.template segment<pow2(NQUBITS - 1)>(pow2(NQUBITS - 1));

            for(int i = 0; i < pow2(NQUBITS - 1); i++) {
                coeffs[i] = low[i] * Gzero + high[i] * Gone;
            }

            return coeffs;
        }
    };

    template<int N, int Q, int NQUBITS, typename T>
    struct SingleGateOpAdder {
        static Vector<pow2(NQUBITS), T> add(const Vector<pow2(NQUBITS), T>& state, const SingleGateCoeffVector<NQUBITS, T>& coeff, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {

        }
    };

    template<int N, int NQUBITS, typename T>
    struct SingleGateOpAdder<N, 0, NQUBITS, T> {
        static Vector<pow2(NQUBITS), T> add(const SingleGateCoeffVector<NQUBITS, T>& coeff, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            return basis<N, NQUBITS - 1>::value() % coeff[N] + SingleGateOpAdder<N - 1, 0, NQUBITS, T>::add(coeff, Gzero, Gone);
        }
    };

    template<int NQUBITS, typename T>
    struct SingleGateOpAdder<0, 0, NQUBITS, T> {
        static Vector<pow2(NQUBITS), T> add(const SingleGateCoeffVector<NQUBITS, T>& coeff, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            return basis<0, NQUBITS - 1>::value() % coeff[0];
        }
    };

    template<int N, int NQUBITS, typename T>
    struct SingleGateOpAdder<N, NQUBITS - 1, NQUBITS, T> {
        static Vector<pow2(NQUBITS), T> add(const SingleGateCoeffVector<NQUBITS, T>& coeff, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            return coeff[N] % basis<N, NQUBITS - 1>::value() + SingleGateOpAdder<N - 1, NQUBITS - 1, NQUBITS, T>::add(coeff, Gzero, Gone);
        }
    };

    template<int NQUBITS, typename T>
    struct SingleGateOpAdder<0, NQUBITS - 1, NQUBITS, T> {
        static Vector<pow2(NQUBITS), T> add(const SingleGateCoeffVector<NQUBITS, T>& coeff, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            return coeff[0] % basis<0, NQUBITS - 1>::value();
        }
    };

    template<int Q, int NQUBITS, typename T>
    struct SingleGateOp {
        static Vector<pow2(NQUBITS), T> calculate(const Vector<pow2(NQUBITS), T>& state, const Eigen::Matrix<Complex<T>, 2, 2>& op) {
            Vector<2, T> Gzero = op * zero<T>();
            Vector<2, T> Gone = op * one<T>();
            SingleGateCoeffVector<NQUBITS, T> coeff = SingleGateCoeffVec<Q, NQUBITS, T>::calculate(state, Gzero, Gone);
            return SingleGateOpAdder<pow2(NQUBITS - 1) - 1, Q, NQUBITS, T>::add(coeff, Gzero, Gone);
        }
    };
};
#endif
