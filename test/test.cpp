#include "../src/Qatux.hpp"
#include <iostream>

using namespace Qatux;
int main(int argc, char **argv) {
    Vector<2, float> a = supPlus<float>();

    Vector<2, float> b = 1.0 * zero<float>() + 2.0 * one<float>();
    Vector<2, float> c = 1.0 * zero<float>() + 3.0 * one<float>();
    State<2> state(a % a);

    state.showOutcomes();

    std::cout << (one<float>() % one<float>()) % one<float>() << std::endl << std::endl;
    std::cout << one<float>() % (one<float>() % one<float>()) << std::endl << std::endl;

    std::cout << c % b << std::endl << std::endl;

    std::cout << zero<float>() % zero<float>() << std::endl << std::endl;
    std::cout << zero<float>() % one<float>() << std::endl << std::endl;
    std::cout << one<float>() % zero<float>() << std::endl << std::endl;
    std::cout << one<float>() % one<float>() << std::endl << std::endl;
    return 0;
}
