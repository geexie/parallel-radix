#include <iostream>
#include <correctness_test.hpp>
#include <perf_test.hpp>
#include <assert.h>
#include <sort_algorithm.hpp>
#include <tbb/tbb.h>
#include <tclap/CmdLine.h>
#include <tclap/ValueArg.h>
#include <tclap/SwitchArg.h>

int main(int argc, char** argv)
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

    res = sortCorrectnessTest<unsigned int, radix::CudaRadix<unsigned int> >(16384, 0, 10000);
    assert(res);
    std::cout << "Radix sort for type: unsigned int passed" << std::endl;

    res = sortCorrectnessTest<int, radix::CudaRadix<int> >(16384, -5000, 5000);
    assert(res);
    std::cout << "Radix sort for type: signed int passed" << std::endl;

    res = sortCorrectnessTest<float, radix::CudaRadix<float> >(64, -256.0f, 256.0f);
    assert(res);
    std::cout << "Radix sort for type: single presition floating point passed" << std::endl;

    res = sortCorrectnessTest<double, radix::CudaRadix<double> >(16384, -256.0, 256.0);
    assert(res);
    std::cout << "Radix sort for type: double presition floating point passed" << std::endl;

    std::cout << "Regression test passed" << std::endl;

    TCLAP::CmdLine cmd("Parallel radix sort", ' ', "0.9");

    TCLAP::SwitchArg perf_cuda("c","cuda","Run CUDA perfomance test", false);
    cmd.add( perf_cuda );

    TCLAP::SwitchArg perf_cilk("i","intel_cilk","Run cilk plus perfomance test", false);
    cmd.add( perf_cilk );

    TCLAP::SwitchArg perf_omp("o","omp","Run OpenMP perfomance test", false);
    cmd.add( perf_omp );

    TCLAP::SwitchArg perf_tbb("t","tbb","Run TBB perfomance test", false);
    cmd.add( perf_tbb );

    TCLAP::SwitchArg perf_sec("s","sec","Run secuental realization perfomance test", false);
    cmd.add( perf_sec );

    TCLAP::ValueArg<size_t> pot_min("m", "min", "Minimal size of tested array  = 2^m", false, 8, "size_t");
    cmd.add( pot_min );

    TCLAP::ValueArg<size_t> pot_max("n","max","Maximal size of tested array  = 2^n",false,26,"size_t");
    cmd.add( pot_max );

    cmd.parse( argc, argv );
    size_t min = pot_min.getValue();
    size_t max = pot_max.getValue();

    std::vector<radix::TestInfo> results;

    if(perf_cuda.getValue())
        results.push_back(sortPerformanceTest<radix::CudaRadix<unsigned>, unsigned>(min, max, 0, 10000));

    if(perf_omp.getValue())
        results.push_back(sortPerformanceTest<radix::OpenMPRadix<unsigned>, unsigned>(min, max, 0, 10000));

    if(perf_tbb.getValue())
        results.push_back(sortPerformanceTest<radix::TBBRadix<unsigned>, unsigned>(min, max, 0, 10000));

    if(perf_sec.getValue())
        results.push_back(sortPerformanceTest<radix::SequentalRadix<unsigned>, unsigned>(min, max, 0, 10000));

    if(perf_cilk.getValue())
        results.push_back(sortPerformanceTest<radix::SilkRadix<unsigned>, unsigned>(min, max, 0, 10000));

    for(size_t i = 0; i < results.size(); i++)
        for(size_t j = 0; j < results[i].measurements.size(); j++)
            std::cout << results[i].measurements[j].first << "\t" << results[i].measurements[j].second << std::endl;

    return 0;
}
