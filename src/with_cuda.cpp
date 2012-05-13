#include <sort_algorithm.hpp>
#include <traits.hpp>
#include <cuda_radix.h>

#include <cuda.h>
#include <cuda_runtime_api.h>

template<typename T>
void radix::CudaRadix<T>::operator ()(T* typed_a, size_t len) const
{
    typedef radix::radix_traits<T>::integer_type I;
    I *a = (I*)typed_a;
    I *d = new I[len];
    I* t;

    unsigned is, id0, id1;

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
    t = a; a = d; d = t;
    delete[] d;
}

template<>
void radix::CudaRadix<unsigned int>::operator ()(unsigned int* a, size_t len) const
{
    unsigned int* d_a, *d_d;
    cudaMalloc((void**)&d_a, len * sizeof(unsigned int));
    cudaMemcpy(d_a, a, len * sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_d, len * sizeof(unsigned int));

    cuda_radix_coller(d_a, d_d, len);

    cudaMemcpy(a, d_a, len * sizeof(unsigned int), cudaMemcpyDeviceToHost);
    cudaFree((void**)d_a);
    cudaFree((void**)&d_d);
}

template void radix::CudaRadix<float>::operator()(float*, size_t) const;
template void radix::CudaRadix<double>::operator()(double*, size_t) const;
template void radix::CudaRadix<signed int>::operator()(signed int*, size_t) const;
