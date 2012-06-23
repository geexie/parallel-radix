#include <sort_algorithm.hpp>
#include <traits.hpp>
#include <cuda_radix.h>

#include <cuda.h>
#include <cuda_runtime_api.h>

template<typename T>
double radix::CudaRadix<T>::operator ()(T* a, size_t len) const
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

template double radix::CudaRadix<float>::operator()(float*, size_t) const;
template double radix::CudaRadix<double>::operator()(double*, size_t) const;
template double radix::CudaRadix<signed int>::operator()(signed int*, size_t) const;
template double radix::CudaRadix<unsigned int>::operator()(unsigned int*, size_t) const;
