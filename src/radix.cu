#include <cuda_radix.h>
#include<iostream>

typedef unsigned int ui;

__global__ void bitonic_sort_step(ui *dev_values, int j, int k)
{
    unsigned int i, ixj; /* Sorting partners: i and ixj */
    i = threadIdx.x + blockDim.x * blockIdx.x;
    ixj = i^j;

    if ((ixj)>i)
    {
        if ((i&k)==0)
        {
            /* Sort ascending */
            if (dev_values[i]>dev_values[ixj])
            {
                /* exchange(i,ixj); */
                float temp = dev_values[i];
                dev_values[i] = dev_values[ixj];
                dev_values[ixj] = temp;
            }
        }
        if ((i&k)!=0)
        {
            /* Sort descending */
            if (dev_values[i]<dev_values[ixj])
            {
                /* exchange(i,ixj); */
                float temp = dev_values[i];
                dev_values[i] = dev_values[ixj];
                dev_values[ixj] = temp;
            }
        }
    }
}

__device__ void blelloch_scan_step(ui *g_idata, ui *g_odata , int n, ui _bit, int shift)
{
    extern __shared__ ui temp1[];
    int thid = threadIdx.x;
    int global_tid = blockIdx.x * blockDim.x * 2;
    ui* bit = &temp1[0];
    ui* idx = &temp1[n];
    size_t tid_offset = 2 * thid;
    {
        int offset = 1;
        bit[2 * thid] = ((g_idata[global_tid + tid_offset] & _bit) >> shift);
        idx[2 * thid] = 1 - bit[tid_offset];
        bit[2 * thid + 1] = ((g_idata[global_tid + tid_offset + 1] & _bit) >> shift);
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

                float t = idx[ai];
                idx[ai] = idx[bi];
                idx[bi] += t;
            }
        }
        __syncthreads();
        ui middle = idx[n - 1] + 1 - bit[n - 1];

        idx[tid_offset]     = (bit[tid_offset + 0])? (tid_offset     - idx[tid_offset + 0] + bit[tid_offset + 0] * middle) : idx[tid_offset + 0];
        idx[tid_offset + 1] = (bit[tid_offset + 1])? (tid_offset + 1 - idx[tid_offset + 1] + bit[tid_offset + 1] * middle) : idx[tid_offset + 1];

        __syncthreads();
        bit[idx[tid_offset + 0]] = g_idata[global_tid + tid_offset + 0];
        bit[idx[tid_offset + 1]] = g_idata[global_tid + tid_offset + 1];
        size_t step  = ((shift == 31) && (blockIdx.x & 1))?blockDim.x * 2 - tid_offset -1: tid_offset;
        size_t step1 = ((shift == 31) && (blockIdx.x & 1))?blockDim.x * 2 - tid_offset -2 : tid_offset + 1;
         __syncthreads();
        g_odata[global_tid + tid_offset + 0] = bit[step];
        g_odata[global_tid + tid_offset + 1] = bit[step1];
    }
}

__global__ void blelloch_scan_radix(ui *a_p, ui *d_p , int n)
{
    int iter = 0;

    for (ui _bit = 1; _bit; _bit <<= 1)
    {
        __syncthreads();
        if(iter % 2 == 0)
            blelloch_scan_step(a_p, d_p, n, _bit, iter++);
        else
            blelloch_scan_step(d_p, a_p, n, _bit, iter++);

    }
}

void cuda_radix_coller(unsigned int* a, unsigned int* d, size_t len)
{
    int threads = 256;
    int log_threads = 8;
    int groups = len >> log_threads;
    blelloch_scan_radix<<<groups, (threads >> 1), 2 * threads * 4>>>(a, d, threads);

    for ( size_t k = threads; k <= len; k <<= 1)
    {
        for (size_t j=k>>1; j>0; j=j>>1)
        {
            bitonic_sort_step<<<groups, threads >>>(a, j, k);
        }
    }
}
