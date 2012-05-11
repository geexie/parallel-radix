#include <cuda_radix.h>

void __global__ radix_kernel(unsigned int* a, unsigned int* d, size_t len)
{
    unsigned int tid = threadIdx.x;
    extern __shared__ unsigned int bounds[];
    unsigned int* = bounds[tid];

}

void cuda_radix_coller(unsigned int* a, unsigned int* d, size_t len)
{
    int threads = 32;
    int groups = 1;
    radix_kernel<groups,threads,(threads * 2)>(a, d, len);
}
