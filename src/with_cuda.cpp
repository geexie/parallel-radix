#include <sort_algorithm.hpp>
#include <traits.hpp>
#include <cuda_radix.h>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>

namespace radix {

template<typename T>
double CudaRadix<T>::operator ()(T* a, size_t len) const
{
    T* d_a, *d_d;
    cudaMalloc((void**)&d_a, len * sizeof(T));
    cudaMemcpy(d_a, a, len * sizeof(T), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_d, len * sizeof(T));

    double time = cuda_radix_coller(d_a, d_d, len);

    cudaMemcpy(a, d_a, len * sizeof(T), cudaMemcpyDeviceToHost);
    cudaFree((void**)d_a);
    cudaFree((void**)&d_d);

    return time;
}

template double CudaRadix<float>::operator()(float*, size_t) const;
template double CudaRadix<double>::operator()(double*, size_t) const;
template double CudaRadix<signed int>::operator()(signed int*, size_t) const;
template double CudaRadix<unsigned int>::operator()(unsigned int*, size_t) const;
}
