#include "../src/Qatux.hpp"
#include <iostream>

using namespace Qatux;
int main(int argc, char **argv) {
    Vector<2, float> a = supPlus<float>();
    State<2> state(a % a);

    state.showOutcomes();

    std::cout << basis<3, 2, float>::value() << std::endl << std::endl;
    return 0;
}
