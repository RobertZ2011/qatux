#include <iostream>
#include "../src/Qatux.hpp"

using namespace Qatux;

int main(int argc, char **argv) {
    State<2> state(Basis<2, 2>::value());

    state.showOutcomes();
    state.hadamard<0>();
    state.cnot<0, 1>();
    state.showOutcomes();
    return 0;
}
