#include <sort_algorithm.hpp>
#include <traits.hpp>
#include <stdio.h>
#include <timing.h>
#include "cilk/cilk.h"

template<typename T>
double radix::SilkRadix<T>::operator ()(T* typed_a, size_t len) const
{
    typedef radix::radix_traits<T>::integer_type I;
    I *a = (I*)typed_a;
    I *d = new I[len];
    I* t;

    size_t block = __cilkrts_get_nworkers();
    size_t* l_bounds = new size_t[block];
    size_t* r_bounds = new size_t[block];

    double time  = (double)radix::getTickCount();
    for (I bit = 1; ( bit & ~radix::radix_traits<T>::MSB_mask); bit <<= 1)
    {
        for(int i = 0; i < block; i++)
        {
            l_bounds[i] = r_bounds[i] = 0;
        }

        cilk_for(int k = 0; k < block; k++)
        {
            for (int i = 0; i < k * len / block; ++i)
            {
                if(!!(((a[i] >> 1) ^ a[i]) & bit))
                    r_bounds[k]++;
                else
                    l_bounds[k]++;
            }
        }

        cilk_for(int k = 0; k < block; k++)
        {
            for (int i = k * len / block; i < (k+1) * len / block; ++i)
            {
                if(!!(((a[i] >> 1) ^ a[i]) & bit))
                    d[len - 1 - r_bounds[k]++] = a[i];
                else
                    d[l_bounds[k]++] = a[i];
            }
        }
        t = a; a = d; d = t;
    }

    //MSD in inverse opder
    size_t l_bound = 0;
    size_t r_bound = 0;

    //counting loop
    for(size_t i = 0; i < len; ++i)
    {
        if (a[i] & radix::radix_traits<T>::MSB_mask) ++r_bound;
    }
    size_t middle = r_bound - 1;

    //permutation loop
    for(size_t i = 0; i < len; ++i)
    {
        if((a[i] & radix::radix_traits<T>::MSB_mask))
        {
            d[radix_traits<T>::index(l_bound, middle)] = a[i];
            l_bound++;
        }
        else
            d[r_bound++] = a[i];
    }
    time = ((double)radix::getTickCount() - time)/radix::getTickFrequency();
    t = a; a = d; d = t;
    delete[] d;
    return time;
}

struct bounds
{
    size_t l_bound;
    size_t r_bound;
};


template<>
double radix::SilkRadix<unsigned int>::operator ()(unsigned int* a, size_t len) const
{
    unsigned int *d = new unsigned int[len];
    unsigned int *t;

    size_t block = __cilkrts_get_nworkers();
    size_t* l_bounds = new size_t[block];
    size_t* r_bounds = new size_t[block];

    double time = (double)radix::getTickCount();
    for (unsigned bit = 1; bit; bit <<= 1)
    {
        for(int i = 0; i < block; i++)
        {
            l_bounds[i] = r_bounds[i] = 0;
        }

        cilk_for(int k = 0; k < block; k++)
        {
            for (int i = 0; i < k * len / block; ++i)
            {
                if(!!(((a[i] >> 1) ^ a[i]) & bit))
                    r_bounds[k]++;
                else
                    l_bounds[k]++;
            }
        }

        cilk_for(int k = 0; k < block; k++)
        {
            for (int i = k * len / block; i < (k+1) * len / block; ++i)
            {
                if(!!(((a[i] >> 1) ^ a[i]) & bit))
                    d[len - 1 - r_bounds[k]++] = a[i];
                else
                    d[l_bounds[k]++] = a[i];
            }
        }
        t = a; a = d; d = t;
    }
    time = ((double)radix::getTickCount() - time)/radix::getTickFrequency();
    delete[] d;

    return time;
}

template double radix::SilkRadix<float>::operator()(float*, size_t) const;
template double radix::SilkRadix<double>::operator()(double*, size_t) const;
template double radix::SilkRadix<signed int>::operator()(signed int*, size_t) const;
