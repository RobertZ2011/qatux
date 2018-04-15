// Let's take a look at a 3 qubit state, written as
// V000|000> + V001|001> + V010|010> + V011|011> +
// V100|100> + V101|101> + V110|110> + V111|111>
// Now apply an operator G to the first qubit
// V000(|00> ⊗ G|0>) + V001(|00> ⊗ G|1>) + V010(|01> ⊗ G|0>) + V011(|01> ⊗ G|1>) +
// V100(|10> ⊗ G|0>) + V101(|10> ⊗ G|1>) + V110(|11> ⊗ G|0>) + V111(|11> ⊗ G|1>)
// No factor out the 2 qubit basis vectors
// |00> ⊗ (V000*G|0> + V001*G|1>) +
// |01> ⊗ (V010*G|0> + V011*G|1>) +
// |10> ⊗ (V100*G|0> + V101*G|1>) +
// |11> ⊗ (V110*G|0> + V111*G|1>)
// A similar method can be used for the second qubit to yield
// |0> ⊗ G|0> ⊗ (V000*|0> + V001*|1>) +
// |0> ⊗ G|1> ⊗ (V010*|0> + V011*|1>) +
// |1> ⊗ G|0> ⊗ (V100*|0> + V101*|1>) +
// |1> ⊗ G|1> ⊗ (V110*|0> + V111*|1>)
// And for the last qubit:
// (V000*G|0> + V100*G|1>) ⊗ |00> +
// (V001*G|0> + V101*G|1>) ⊗ |01> +
// (V010*G|0> + V110*G|1>) ⊗ |10> +
// (V011*G|0> + V111*G|1>) ⊗ |11>
// This is how we'll perform single qubit gates

#ifndef QATUX_SINGLE_GATE
#define QATUX_SINGLE_GATE

namespace Qatux {
    template<int NQUBITS, typename T>
    using SingleGateCoeffVector = Eigen::Array<Vector<2, T>, pow2(NQUBITS - 1), 1>;

    // Coefficient vector was the best name I could come up with for the vector
    // that's calculated using the state coefficients Vx
    // this is the general version
    // Gzero and Gone are G|0> and G|1>
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

    // specialized template for qubit 0
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

    // specialized template for qubit NQUBITS - 1
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

    // Used to add together the vectors to produce the new state
    // General template for Q != 0, Q != 1 and Q != NQUBITS - 1
    template<int N, int Q, int NQUBITS, typename T>
    struct SingleGateOpAdder {
        static Vector<pow2(NQUBITS), T> add(const SingleGateCoeffVector<NQUBITS, T>& coeff, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            constexpr int upperBasis = N >> (Q + 1);     //gets the bits above bit Q
            constexpr int upperSize = NQUBITS - (Q + 1); //number of bits in upperBasis
            constexpr int lowerBasis = ~(1 << Q) & N; //gets the bits below bit Q
            constexpr int lowerSize = Q;              //number of bits in lowerBasis
            auto selectedQubit = (N & 0x1) ? Gone : Gzero;
            return basis<upperBasis, upperSize>::value () % selectedQubit % basis<lowerBasis, lowerSize>::value() % coeff[N] + SingleGateOpAdder<N - 1, Q, NQUBITS, T>::add(coeff, Gzero, Gone);
        }
    };

    // specialized template to end recursion for Q != 0, Q != 1 and Q != NQUBITS - 1
    template<int Q, int NQUBITS, typename T>
    struct SingleGateOpAdder<0, Q, NQUBITS, T> {
        static Vector<pow2(NQUBITS), T> add(const SingleGateCoeffVector<NQUBITS, T>& coeff, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            constexpr int upperSize = NQUBITS - (Q + 1); //number of bits in upperBasis
            constexpr int lowerSize = Q;                 //number of bits in lowerBasis
            return basis<0, upperSize>::value () % Gzero % basis<0, lowerSize>::value() % coeff[0];
        }
    };

    // specialized template for qubit 0
    template<int N, int NQUBITS, typename T>
    struct SingleGateOpAdder<N, 0, NQUBITS, T> {
        static Vector<pow2(NQUBITS), T> add(const SingleGateCoeffVector<NQUBITS, T>& coeff, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            return basis<N, NQUBITS - 1>::value() % coeff[N] + SingleGateOpAdder<N - 1, 0, NQUBITS, T>::add(coeff, Gzero, Gone);
        }
    };

    // specialized template for qubit 0 that ends the recursion
    template<int NQUBITS, typename T>
    struct SingleGateOpAdder<0, 0, NQUBITS, T> {
        static Vector<pow2(NQUBITS), T> add(const SingleGateCoeffVector<NQUBITS, T>& coeff, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            return basis<0, NQUBITS - 1>::value() % coeff[0];
        }
    };

    // specialized template for qubit 1
    template<int N, int NQUBITS, typename T>
    struct SingleGateOpAdder<N, 1, NQUBITS, T> {
        static Vector<pow2(NQUBITS), T> add(const SingleGateCoeffVector<NQUBITS, T>& coeff, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            constexpr int upperBasis = N >> 2;     //gets the bits above bit Q
            constexpr int upperSize = NQUBITS - 2; //number of bits in upperBasis
            auto selectedQubit = (N & 0x1) ? Gone : Gzero;
            return basis<upperBasis, upperSize>::value () % selectedQubit % coeff[N] + SingleGateOpAdder<N - 1, 1, NQUBITS, T>::add(coeff, Gzero, Gone);
        }
    };

    // specialized template for qubit 1 that ends the recursion
    template<int NQUBITS, typename T>
    struct SingleGateOpAdder<0, 1, NQUBITS, T> {
        static Vector<pow2(NQUBITS), T> add(const SingleGateCoeffVector<NQUBITS, T>& coeff, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            constexpr int upperSize = NQUBITS - 2; //number of bits in upperBasis
            return basis<0, upperSize>::value () % Gzero % coeff[0];
        }
    };

    // specialized template for NQUBITS = 2, Q = 1
    template<int N, typename T>
    struct SingleGateOpAdder<N, 1, 2, T> {
        static Vector<4, T> add(const SingleGateCoeffVector<2, T>& coeff, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            return coeff[N] % basis<N, 1>::value() + SingleGateOpAdder<N - 1, 1, 2, T>::add(coeff, Gzero, Gone);
        }
    };

    // specialized template for NQUBITS = 2, Q = 1 that ends the recursion
    template<typename T>
    struct SingleGateOpAdder<0, 1, 2, T> {
        static Vector<4, T> add(const SingleGateCoeffVector<2, T>& coeff, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            return coeff[0] % basis<0, 1>::value();
        }
    };

    // specialized template for qubit NQUBITS - 1
    template<int N, int NQUBITS, typename T>
    struct SingleGateOpAdder<N, NQUBITS - 1, NQUBITS, T> {
        static Vector<pow2(NQUBITS), T> add(const SingleGateCoeffVector<NQUBITS, T>& coeff, const Vector<2, T>& Gzero, const Vector<2, T>& Gone) {
            return coeff[N] % basis<N, NQUBITS - 1>::value() + SingleGateOpAdder<N - 1, NQUBITS - 1, NQUBITS, T>::add(coeff, Gzero, Gone);
        }
    };

    // specialized template for qubit NQUBITS - 1 that ends the recursion
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
