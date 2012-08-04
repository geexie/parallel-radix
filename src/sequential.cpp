#include <sort_algorithm.hpp>
#include <traits.hpp>
#include <timing.h>

namespace radix
{

template<typename T>
double SequentalRadix<T>::operator ()(T* typed_a, size_t len) const
{
    typedef typename radix::radix_traits<T>::integer_type I;
    I *a = (I*)typed_a;
    I *d = new I[len];
    I* t;

    unsigned is, id0, id1;
    double time = (double)radix::getTickCount();
    for (I bit = 1; (bit & ~radix::radix_traits<T>::MSB_mask); bit <<= 1, t = a, a = d, d = t)
    {
        for (is = id0 = 0, id1 = len; id1 > id0; ++is)
            //reopder by Gray order
            if (((a[is] >> 1) ^ a[is]) & bit)
                d[--id1] = a[is];
            else
                d[id0++] = a[is];
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

template<>
double SequentalRadix<unsigned int>::operator ()(unsigned int* a, size_t len) const
{
        unsigned int *s = a;
        unsigned int *d = new unsigned int[len];
        unsigned int *t;
        unsigned bit, is, id0, id1;
        double time = (double)radix::getTickCount();
        for (bit = 1; bit; bit <<= 1, t = s, s = d, d = t)
            for (is = id0 = 0, id1 = len; id1 > id0; ++is)
                //reopder by Gray order
                if (((s[is] >> 1) ^ s[is]) & bit)
                    d[--id1] = s[is];
                else
                    d[id0++] = s[is];
        time = ((double)radix::getTickCount() - time)/radix::getTickFrequency();
        delete[] d;

        return time;
}

template double SequentalRadix<float>::operator()(float*, size_t) const;
template double SequentalRadix<double>::operator()(double*, size_t) const;
template double SequentalRadix<signed int>::operator()(signed int*, size_t) const;

}
