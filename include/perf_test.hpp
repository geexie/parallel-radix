#ifndef __perf_test_hpp__
#define __perf_test_hpp__

#include <traits.hpp>

template<typename Alg, typename T>
std::pair<size_t, double> sortPerformanceTestStep(size_t size, T l_bound, T r_bound)
{
    T * arr =  new T[size];

    T sign = random_traits<T>::sign;
    random_traits<T>::urand_gen urand(l_bound, r_bound);
    std::tr1::mt19937 eng;

    for (size_t i = 0; i < size; ++i)
    {
        *(arr + i) = sign * urand(eng);
        sign *= random_traits<T>::sign;
    }

    Alg alg;
    return std::make_pair(size, alg(arr, size));
}

template<typename Alg, typename T>
radix::TestInfo sortPerformanceTest(size_t pot_min, size_t pot_max,T l_bound, T r_bound)
{
    radix::TestInfo res;
    for(size_t i = 1U << pot_min; i <= 1U << pot_max; i<<=1)
    {
        res.add(sortPerformanceTestStep<Alg, T>(i, l_bound, r_bound));
    }
    return res;
}

#endif
