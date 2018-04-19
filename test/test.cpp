#include <iostream>
#include "../src/Qatux.hpp"

using namespace Qatux;

int main(int argc, char **argv) {
    State<float> state(4);

    state = basis(0, 4);

    state.showOutcomes();
    for(int i = 0; i < 10000; i++) {
        state.hadamard(0);
        state.hadamard(1);
        state.hadamard(2);
        state.hadamard(3);
    }
    state.showOutcomes();
    return 0;
}
