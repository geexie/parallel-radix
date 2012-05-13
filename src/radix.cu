#include <cuda_radix.h>
#include<iostream>

__device__ void blelloch_scan_step(unsigned int *g_idata, unsigned int *g_odata , int n, unsigned int _bit, int shift)
{
    extern __shared__ unsigned int temp1[];
    int thid = threadIdx.x;
    unsigned int* bit = &temp1[0];
    unsigned int* idx = &temp1[n];
    size_t tid_offset = 2 * thid;

    {
        int offset = 1;
        bit[2 * thid] = (g_idata[tid_offset] & _bit) >> shift;
        idx[2 * thid] = 1 - bit[tid_offset];
        bit[2 * thid + 1] = (g_idata[tid_offset + 1] & _bit) >> shift;
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
        unsigned int middle = idx[n - 1] + 1 - bit[n - 1];

        idx[tid_offset]     = (bit[tid_offset + 0])? (tid_offset     - idx[tid_offset + 0] + bit[tid_offset + 0] * middle) : idx[tid_offset + 0];
        idx[tid_offset + 1] = (bit[tid_offset + 1])? (tid_offset + 1 - idx[tid_offset + 1] + bit[tid_offset + 1] * middle) : idx[tid_offset + 1];

        g_odata[idx[tid_offset]]     = g_idata[tid_offset];
        g_odata[idx[tid_offset + 1]] = g_idata[tid_offset + 1];
    }
}

__global__ void blelloch_scan_radix(unsigned int *a_p, unsigned int *d_p , int n)
{
    int iter = 0;

    for (unsigned int _bit = 1; _bit; _bit <<= 1)
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
    int threads = len >> 1;
    int groups = 1;
    blelloch_scan_radix<<<groups, threads, 2 * len * 4>>>(a, d, len);
}
