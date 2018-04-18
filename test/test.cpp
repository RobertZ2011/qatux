#include <iostream>
#include "../src/Qatux.hpp"

using namespace Qatux;

int main(int argc, char **argv) {
    State<7> state(Basis<0, 7>::value());

    state.showOutcomes();
    state.hadamard<0>();
    state.hadamard<1>();
    state.hadamard<2>();
    state.hadamard<3>();
    state.showOutcomes();
    return 0;
}
