#include <cuda_radix.h>
#include <iostream>
#include <traits.hpp>
#include <timer.cuh>
#include <stdio.h>

// ========================= botonic butterfly function ======================== //
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

// ========================= work efficient blelloch sort ====================== //
template<typename T, typename I>
__device__ void blelloch_scan_step(T *g_idata, T *g_odata , size_t n, I _bit, int shift, volatile unsigned int* bit, unsigned int* idx)
{
    int thid = threadIdx.x;
    int global_tid = blockIdx.x * blockDim.x * 2;

    size_t tid_offset = 2 * thid;

    int offset = 1;

    bit[2 * thid + 0] = ((I)(radix::radix_traits<T>::as_integer(g_idata[global_tid + tid_offset]) & _bit) >> shift);
    bit[2 * thid + 1] = ((I)(radix::radix_traits<T>::as_integer(g_idata[global_tid + tid_offset + 1]) & _bit) >> shift);

    idx[tid_offset + 0] = 1 - bit[tid_offset + 0];
    idx[tid_offset + 1] = 1 - bit[tid_offset + 1];

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

    if(shift == radix::radix_traits<T>::MSB_bit)
    {
        g_odata[global_tid + radix::radix_traits<T>::index(idx[tid_offset + 0], middle)] = g_idata[global_tid + tid_offset + 0];
        g_odata[global_tid + radix::radix_traits<T>::index(idx[tid_offset + 1], middle)] = g_idata[global_tid + tid_offset + 1];
    }
    else
    {
        g_odata[global_tid + idx[tid_offset + 0]] = g_idata[global_tid + tid_offset + 0];
        g_odata[global_tid + idx[tid_offset + 1]] = g_idata[global_tid + tid_offset + 1];
    }
}

template<typename T>
__global__ void blelloch_scan_radix(T *a_p, T *d_p , int len)
{
    int shift = 0;
    typedef typename radix::radix_traits<T>::integer_type I;

    extern __shared__ char temp1[];

    unsigned int* bit = (unsigned int*)temp1;
    unsigned int* idx = (unsigned int*)(temp1 + len * sizeof(unsigned int));

    for (I _bit = 1; _bit; _bit <<= 1)
    {
        __syncthreads();
        if(shift % 2 == 0)
            blelloch_scan_step<T, I>(a_p, d_p, len, _bit, shift++, bit, idx);
        else
            blelloch_scan_step<T, I>(d_p, a_p, len, _bit, shift++, bit, idx);
    }
}
// ========================= ballot streem compression based radix ============= //
__device__ __forceinline__ unsigned int __bfi(unsigned int x, unsigned int y, unsigned int bit, unsigned int numBits)
{
    unsigned int ret;
    asm("bfi.b32 %0, %1, %2, %3, %4;" :
    "=r"(ret) : "r"(y), "r"(x), "r"(bit), "r"(numBits));
    return ret;
}

template<typename I>
__device__ __forceinline__ I grey(I x)
{
    return (x >> 1) ^ x;
}

// typedef unsigned int uint;

template<typename T, typename I, int num_warps, int log_warps>
__device__ void ballot_scan(const T* a, T* d, I _bit, int shift)
{
    unsigned int tid = threadIdx.x;
    unsigned int lane = 31 & tid;
    unsigned int warp = tid / 32;
    T ini_val = a[blockDim.x * blockIdx.x + tid];

    // ones
    T val = ini_val;
    uint flag = ( grey(radix::radix_traits<T>::as_integer(val))  & _bit) >> shift;

    unsigned int bits = __ballot(flag);
    unsigned int mask = __bfi(0, 0xffffffff, 0, lane);
    unsigned int exc  = __popc(mask & bits);
    unsigned int warpTotal = __popc(bits);

    __shared__ volatile unsigned int totals1[num_warps];
    __shared__ volatile unsigned int totals2[num_warps];
    if (!lane) totals1[warp] = warpTotal;

    // Inclusive scan the warp totals.
    __syncthreads();

    if(tid < num_warps)
    {
        unsigned int x = totals1[tid];
        for(int i = 0; i < log_warps; ++i)
        {
            uint offset = 1<< i;
            if (tid >= offset)
                x += totals1[tid - offset];
            totals1[tid] = x;
        }
    }
    __syncthreads();

    // // Add the scanned warp totals into exc.
    unsigned int blockTotal = totals1[num_warps - 1];
    exc += totals1[warp] - warpTotal;

    __shared__ volatile T shared1[num_warps * 32];
    __shared__ volatile T shared2[num_warps * 32];

    if(flag) shared1[exc] = val;
    __syncthreads();

    if(tid < blockTotal)
    {
        val = shared1[tid];
        d[blockDim.x * blockIdx.x + (blockDim.x - 1 - tid)] = val;
    }

    __syncthreads();
    // zerros
    val = ini_val;
    flag = (1 - (( grey(radix::radix_traits<T>::as_integer(val)) & _bit) >> shift));

    bits = __ballot(flag);
    mask = __bfi(0, 0xffffffff, 0, lane);
    exc  = __popc(mask & bits);
    warpTotal = __popc(bits);

    if (!lane) totals2[warp] = warpTotal;

    // Inclusive scan the warp totals.
    __syncthreads();

    if(tid < num_warps)
    {
        unsigned int x = totals2[tid];
        for(int i = 0; i < log_warps; ++i)
        {
            uint offset = 1<< i;
            if (tid >= offset)
                x += totals2[tid - offset];
            totals2[tid] = x;
        }
    }
    __syncthreads();

    // Add the scanned warp totals into exc.
    blockTotal = totals2[num_warps - 1];
    exc += totals2[warp] - warpTotal;

    if (flag) shared2[exc] = val;
    __syncthreads();

    if(tid < blockTotal)
    {
        val = shared2[tid];
        d[blockDim.x * blockIdx.x + tid] = val;
    }
}

template<typename T, typename I, int num_warps, int log_warps>
__device__ void ballot_scan_des(const T* a, T* d, I _bit, int shift)
{
    unsigned int tid = threadIdx.x;
    unsigned int lane = 31 & tid;
    unsigned int warp = tid / 32;
    T ini_val = a[blockDim.x * blockIdx.x + tid];

    // ones
    T val = ini_val;
    uint flag = 1 - (( grey(radix::radix_traits<T>::as_integer(val))  & _bit) >> shift);

    unsigned int bits = __ballot(flag);
    unsigned int mask = __bfi(0, 0xffffffff, 0, lane);
    unsigned int exc  = __popc(mask & bits);
    unsigned int warpTotal = __popc(bits);

    __shared__ volatile unsigned int totals1[num_warps];
    __shared__ volatile unsigned int totals2[num_warps];
    if (!lane) totals1[warp] = warpTotal;

    // Inclusive scan the warp totals.
    __syncthreads();

    if(tid < num_warps)
    {
        unsigned int x = totals1[tid];
        for(int i = 0; i < log_warps; ++i)
        {
            uint offset = 1<< i;
            if (tid >= offset)
                x += totals1[tid - offset];
            totals1[tid] = x;
        }
    }
    __syncthreads();

    // // Add the scanned warp totals into exc.
    unsigned int blockTotal = totals1[num_warps - 1];
    exc += totals1[warp] - warpTotal;

    __shared__ volatile T shared1[num_warps * 32];
    __shared__ volatile T shared2[num_warps * 32];

    if(flag) shared1[exc] = val;
    __syncthreads();

    if(tid < blockTotal)
    {
        val = shared1[blockTotal - 1 - tid];
        d[blockDim.x * blockIdx.x + (blockDim.x - 1 - tid)] = val;
    }

    __syncthreads();
    // zerros
    val = ini_val;
    flag = ((( grey(radix::radix_traits<T>::as_integer(val)) & _bit) >> shift));

    bits = __ballot(flag);
    mask = __bfi(0, 0xffffffff, 0, lane);
    exc  = __popc(mask & bits);
    warpTotal = __popc(bits);

    if (!lane) totals2[warp] = warpTotal;

    // Inclusive scan the warp totals.
    __syncthreads();

    if(tid < num_warps)
    {
        unsigned int x = totals2[tid];
        for(int i = 0; i < log_warps; ++i)
        {
            uint offset = 1<< i;
            if (tid >= offset)
                x += totals2[tid - offset];
            totals2[tid] = x;
        }
    }
    __syncthreads();

    // Add the scanned warp totals into exc.
    blockTotal = totals2[num_warps - 1];
    exc += totals2[warp] - warpTotal;

    if (flag) shared2[exc] = val;
    __syncthreads();

    if(tid < blockTotal)
    {
        val = shared2[blockTotal - 1 - tid];
        d[blockDim.x * blockIdx.x + tid] = val;
    }
}

template<typename T, int num_warps, int log_warps>
__global__ void ballot_scan_radix(T *a_p, T *d_p , int len)
{
    int shift = 0;
    typedef typename radix::radix_traits<T>::integer_type I;

    if (blockIdx.x & 1)
    {
        for (I _bit = 1; _bit; _bit <<= 1)
        {
            __syncthreads();
            if(shift % 2 == 0)
                ballot_scan_des<T, I, num_warps, log_warps>(a_p, d_p, _bit, shift++);
            else
                ballot_scan_des<T, I, num_warps, log_warps>(d_p, a_p, _bit, shift++);
            // if (threadIdx.x == 0) printf("a - step %u\n", _bit);
        }
    }
    else
    {

        for (I _bit = 1; _bit; _bit <<= 1)
        {
            __syncthreads();
            if(shift % 2 == 0)
                ballot_scan<T, I, num_warps, log_warps>(a_p, d_p, _bit, shift++);
            else
                ballot_scan<T, I, num_warps, log_warps>(d_p, a_p, _bit, shift++);
            // if (threadIdx.x == 0) printf("d - step %u\n", _bit);
        }
    }
}

template<typename T>
double cuda_radix_coller_internal(T* a, T* d, size_t len)
{
    int threads = (len <= 512) ? len : 512;
    int log_threads = 9;
    int groups = (len <= 512) ? 1 : len >> log_threads;
    Timer timer;
    printf("threds %d, groups %d\n", threads, groups);
    timer.go();
    // blelloch_scan_radix<T><<<groups, (threads >> 1), 2 * threads * sizeof(unsigned int)>>>(a, d, threads);
    switch (len)
    {
        case 32  : ballot_scan_radix<T, 1,  0><<<1,      32>>>(a, d, threads);  break;
        case 64  : ballot_scan_radix<T, 2,  1><<<1,      64>>>(a, d, threads);  break;
        case 128 : ballot_scan_radix<T, 4,  2><<<1,      128>>>(a, d, threads); break;
        case 256 : ballot_scan_radix<T, 8,  3><<<1,      256>>>(a, d, threads); break;
        case 512 : ballot_scan_radix<T, 16, 4><<<groups, 512>>>(a, d, threads); break;
        default  : ballot_scan_radix<T, 16, 4><<<groups, 512>>>(a, d, threads); break;
    }
    cudaDeviceSynchronize();
    float ms = timer.measure();
    if(len > 512)
        for (size_t k = threads; k <= len; k <<= 1)
        {
            for (size_t j = k>>1; j > 0; j = j>>1)
            {
                bitonic_sort_step<T><<<groups, threads >>>(a, j, k);
            }
        }
    return ms ;/// 1000.0;
}

double cuda_radix_coller(unsigned int* a, unsigned int* d, size_t len)
{
    printf("HERE\n");
    return cuda_radix_coller_internal<unsigned int>(a,d,len);
}

double cuda_radix_coller(signed int* a, signed int* d, size_t len)
{
    //return cuda_radix_coller_internal<signed int>(a,d,len);
    return 0.0;
}

double cuda_radix_coller(float* a, float* d, size_t len)
{
    return cuda_radix_coller_internal<float>(a,d,len);
}

double cuda_radix_coller(double* a, double* d, size_t len)
{
    //return cuda_radix_coller_internal<double>(a,d,len);
    return 0.0;
}
