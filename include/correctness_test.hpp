#ifndef __correctness_test_hpp__
#define __correctness_test_hpp__

#ifdef _MSC_VER
# include <random>
#else
# include <tr1/random>
#endif

#include <sort_algorithm.hpp>
#include <bitset>

template<typename T>
struct random_traits{};

template<>
struct random_traits<float>
{
    static const signed int sign = -1;
    typedef std::tr1::uniform_real<float> urand_gen;
};

template<>
struct random_traits<double>
{
    static const signed int sign = -1;
    typedef std::tr1::uniform_real<double> urand_gen;
};

template<>
struct random_traits<signed int>
{
    static const signed int sign = -1;
    typedef std::tr1::uniform_int<signed int> urand_gen;
};

template<>
struct random_traits<unsigned int>
{
    static const unsigned int sign = 1;
    typedef std::tr1::uniform_int<unsigned int> urand_gen;
};

template<typename T, typename Alg>
bool sortCorrectnessTest(size_t size, T l_bound, T r_bound)
{
    T * arr =  new T[size];
    T * gold = new T[size];

    T sign = random_traits<T>::sign;
    typename random_traits<T>::urand_gen urand(l_bound, r_bound);
    std::tr1::mt19937 eng;

    for (size_t i = 0; i < size; ++i)
    {
        *(arr + i) = sign * urand(eng);
        *(gold + i) = *(arr + i);
        sign *= random_traits<T>::sign;
    }

    radix::SequentalSTD<T> stds;
    stds(gold, size);

    Alg alg;
    alg(arr, size);
    bool result = true;

    for (size_t i =0; i < size; ++i)
    {
        if(*(gold + i) != *(arr + i))
        {
            std::cout <<"\t" << std::bitset<32>(((*(int*)(gold + i)) >> 1) ^ (*(int*)(gold + i)))
                      << "\t" "expected " << (*(gold + i))
                      <<"\t" << std::bitset<32>(((*(int*)(arr + i)) >> 1) ^ (*(int*)(arr + i)))
                      << "\t" <<"\tactual " << (*(arr + i))<< "\tfor " << i << std::endl;
            result = false;
        }
    }

    return result;
}
#endif
