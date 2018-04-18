#include <iostream>
#include "../src/Qatux.hpp"

using namespace Qatux;

int main(int argc, char **argv) {
    State<4> state(basis<15, 4>::value());

    for(int i = 0; i < 100000; i ++) {
        state.hadamard<2>();
    }
    return 0;
}
