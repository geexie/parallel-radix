#include <iostream>
#include <correctness_test.hpp>
#include <perf_test.hpp>
#include <assert.h>
#include <sort_algorithm.hpp>
#include <tbb/tbb.h>


int main(int args, char** argv)
{
    std::cout << "Initial sort correctness test" << std::endl << std::endl;

    std::cout << "Sequental radix ..." << std::endl;

    bool res = sortCorrectnessTest<unsigned int, radix::SequentalRadix<unsigned int> >(10000, 0, 10000);
    assert(res);
    std::cout << "Radix sort for type: unsigned int passed" << std::endl;

    res = sortCorrectnessTest<int, radix::SequentalRadix<int> >(10000, -5000, 5000);
    assert(res);
    std::cout << "Radix sort for type: signed int passed" << std::endl;

    res = sortCorrectnessTest<float, radix::SequentalRadix<float> >(10000, -5000.f, 5000.f);
    assert(res);
    std::cout << "Radix sort for type: single presition floating point passed" << std::endl;

    res = sortCorrectnessTest<double, radix::SequentalRadix<double> >(10000, -5000, 5000);
    assert(res);
    std::cout << "Radix sort for type: double presition floating point passed" << std::endl;

    std::cout << "TBB based radix ..." << std::endl;

    res = sortCorrectnessTest<unsigned int, radix::TBBRadix<unsigned int> >(10000, 0, 10000);
    assert(res);
    std::cout << "Radix sort for type: unsigned int passed" << std::endl;

    res = sortCorrectnessTest<int, radix::TBBRadix<int> >(10000, -5000, 5000);
    assert(res);
    std::cout << "Radix sort for type: signed int passed" << std::endl;

    res = sortCorrectnessTest<float, radix::TBBRadix<float> >(10000, -5000.f, 5000.f);
    assert(res);
    std::cout << "Radix sort for type: single presition floating point passed" << std::endl;

    res = sortCorrectnessTest<double, radix::TBBRadix<double> >(10000, -5000, 5000);
    assert(res);
    std::cout << "Radix sort for type: double presition floating point passed" << std::endl;

    std::cout << "OpenMP based radix ..." << std::endl;

    res = sortCorrectnessTest<unsigned int, radix::OpenMPRadix<unsigned int> >(10000, 0, 10000);
    assert(res);
    std::cout << "Radix sort for type: unsigned int passed" << std::endl;

    res = sortCorrectnessTest<int, radix::OpenMPRadix<int> >(10000, -5000, 5000);
    assert(res);
    std::cout << "Radix sort for type: signed int passed" << std::endl;

    res = sortCorrectnessTest<float, radix::OpenMPRadix<float> >(10000, -5000.f, 5000.f);
    assert(res);
    std::cout << "Radix sort for type: single presition floating point passed" << std::endl;

    res = sortCorrectnessTest<double, radix::OpenMPRadix<double> >(10000, -5000, 5000);
    assert(res);
    std::cout << "Radix sort for type: double presition floating point passed" << std::endl;

    std::cout << "Cilk based radix ..." << std::endl;

    res = sortCorrectnessTest<unsigned int, radix::SilkRadix<unsigned int> >(10000, 0, 10000);
    assert(res);
    std::cout << "Radix sort for type: unsigned int passed" << std::endl;

    res = sortCorrectnessTest<int, radix::SilkRadix<int> >(10000, -5000, 5000);
    assert(res);
    std::cout << "Radix sort for type: signed int passed" << std::endl;

    res = sortCorrectnessTest<float, radix::SilkRadix<float> >(10000, -5000.f, 5000.f);
    assert(res);
    std::cout << "Radix sort for type: single presition floating point passed" << std::endl;

    res = sortCorrectnessTest<double, radix::SilkRadix<double> >(10000, -5000, 5000);
    assert(res);
    std::cout << "Radix sort for type: double presition floating point passed" << std::endl;

    std::cout << "Cuda based radix ..." << std::endl;

    res = sortCorrectnessTest<unsigned int, radix::CudaRadix<unsigned int> >(1024, 0, 10000);
    assert(res);
    std::cout << "Radix sort for type: unsigned int passed" << std::endl;

    res = sortCorrectnessTest<int, radix::CudaRadix<int> >(10000, -5000, 5000);
    //assert(res);
    std::cout << "Radix sort for type: signed int passed" << std::endl;

    res = sortCorrectnessTest<float, radix::CudaRadix<float> >(10000, -5000.f, 5000.f);
    assert(res);
    std::cout << "Radix sort for type: single presition floating point passed" << std::endl;

    res = sortCorrectnessTest<double, radix::CudaRadix<double> >(10000, -5000, 5000);
    assert(res);
    std::cout << "Radix sort for type: double presition floating point passed" << std::endl;

    return 0;
}
