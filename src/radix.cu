#include <cuda_radix.h>
#include <iostream>
#include <traits.hpp>


template<typename T>
__global__ void bitonic_sort_step(T *a, int j, int k)
{
    size_t i = threadIdx.x + blockDim.x * blockIdx.x;
    size_t ixj = i^j;

    if (ixj > i)
    {
        if ((i & k) == 0)
        {
            if (a[i] > a[ixj])
            {
                T temp = a[i];
                a[i] = a[ixj];
                a[ixj] = temp;
            }
        }
        if ((i & k) != 0)
        {
            if (a[i] < a[ixj])
            {
                T temp = a[i];
                a[i] = a[ixj];
                a[ixj] = temp;
            }
        }
    }
}

template<typename T>
__device__ void blelloch_scan_step(T *g_idata, T *g_odata , size_t n,typename radix::radix_traits<T>::integer_type _bit, int shift)
{
    typedef typename radix::radix_traits<T>::integer_type I;
    extern __shared__ char temp1[];

    int thid = threadIdx.x;
    int global_tid = blockIdx.x * blockDim.x * 2;

    I* bit = &((I*)temp1)[0];
    I* idx = &((I*)temp1)[n];
    size_t tid_offset = 2 * thid;

    int offset = 1;

    bit[2 * thid + 0] = ((radix::radix_traits<T>::as_integer(g_idata[global_tid + tid_offset]) & _bit) >> shift);
    bit[2 * thid + 1] = ((radix::radix_traits<T>::as_integer(g_idata[global_tid + tid_offset + 1]) & _bit) >> shift);

    idx[2 * thid + 0] = 1 - bit[tid_offset + 0];
    idx[2 * thid + 1] = 1 - bit[tid_offset + 1];

    __syncthreads();
    for (int d = (n >> 1); d > 0; d >>= 1)
    {
        __syncthreads();
        if (thid < d)
        {
            int ai = tid_offset * offset + (offset - 1);
            int bi = tid_offset * offset + (offset - 1) + offset;
            idx[bi] += idx[ai];
        }

        offset *= 2;
    }

    __syncthreads();
    if (thid == 0) { idx[n - 1] = 0; }

    for (int d = 1; d < n; d *= 2)
    {
        offset >>= 1;
        __syncthreads();
        if (thid < d)
        {
            int ai = offset * tid_offset + 1 * offset - 1;
            int bi = offset * tid_offset + 2 * offset - 1;

            T t = idx[ai];
            idx[ai] = idx[bi];
            idx[bi] += t;
        }
    }
    __syncthreads();
    int middle = idx[n - 1] + 1 - bit[n - 1];

    idx[tid_offset + 0] = (bit[tid_offset + 0])? (tid_offset + 0 - idx[tid_offset + 0] + middle) : idx[tid_offset + 0];
    idx[tid_offset + 1] = (bit[tid_offset + 1])? (tid_offset + 1 - idx[tid_offset + 1] + middle) : idx[tid_offset + 1];

    __syncthreads();
    T* buffer = &((T*)temp1)[0];
    buffer[idx[tid_offset + 0]] = g_idata[global_tid + tid_offset + 0];
    buffer[idx[tid_offset + 1]] = g_idata[global_tid + tid_offset + 1];

    __syncthreads();

    if(shift == radix::radix_traits<T>::MSB_bit)
    {
        g_odata[global_tid + radix::radix_traits<T>::index(tid_offset + 0, middle)] = buffer[tid_offset + 0];
        g_odata[global_tid + radix::radix_traits<T>::index(tid_offset + 1, middle)] = buffer[tid_offset + 1];
    }
    else
    {
        g_odata[global_tid + tid_offset + 0] = buffer[(tid_offset + 0)];
        g_odata[global_tid + tid_offset + 1] = buffer[(tid_offset + 1)];
    }
}

template<typename T>
__global__ void blelloch_scan_radix(T *a_p, T *d_p , int len)
{
    int shift = 0;
    typedef typename radix::radix_traits<T>::integer_type I;
    for (I _bit = 1; _bit; _bit <<= 1)
    {
        __syncthreads();
        if(shift % 2 == 0)
            blelloch_scan_step<T>(a_p, d_p, len, _bit, shift++);
        else
            blelloch_scan_step<T>(d_p, a_p, len, _bit, shift++);

    }
}

template<typename T>
void cuda_radix_coller_internal(T* a, T* d, size_t len)
{
    int threads = (len <= 256) ? len : 256;
    int log_threads = 8;
    int groups = (len <= 256) ? 1 : len >> log_threads;
    blelloch_scan_radix<T><<<groups, (threads >> 1), 2 * threads * sizeof(T)>>>(a, d, threads);

    if(len > 256)
        for (size_t k = threads; k <= len; k <<= 1)
        {
            for (size_t j = k>>1; j > 0; j = j>>1)
            {
                bitonic_sort_step<T><<<groups, threads >>>(a, j, k);
            }
        }
}

void cuda_radix_coller(unsigned int* a, unsigned int* d, size_t len)
{
    cuda_radix_coller_internal<unsigned int>(a,d,len);
}

void cuda_radix_coller(signed int* a, signed int* d, size_t len)
{
    cuda_radix_coller_internal<signed int>(a,d,len);
}

void cuda_radix_coller(float* a, float* d, size_t len)
{
    cuda_radix_coller_internal<float>(a,d,len);
}

void cuda_radix_coller(double* a, double* d, size_t len)
{
    cuda_radix_coller_internal<double>(a,d,len);
}
