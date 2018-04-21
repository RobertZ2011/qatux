#include <iostream>
#include "../src/Qatux.hpp"

using namespace Qatux;

int main(int argc, char **argv) {
    State<float> state(5, basis(8, 5));

    state.showOutcomes();

    state.cnot(1, 3);
    state.showOutcomes();
    return 0;
}
