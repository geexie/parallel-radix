#include <sort_algorithm.hpp>
#include <traits.hpp>
#include <omp.h>

#include <iostream>

template<typename T>
radix::OpenMPRadix<T>::OpenMPRadix()
{
    size_t num = 0;
#pragma omp parallel default(none) shared(num)
    {
        num = omp_get_num_threads();
    }

    num_threads = num;
    l_partial =   new unsigned int[num_threads];
    r_partial =  new unsigned int[num_threads];
    l_sum =      new unsigned int[num_threads];
    r_sum =     new unsigned int[num_threads];
}

template<typename T>
radix::OpenMPRadix<T>::~OpenMPRadix()
{
    delete[] l_partial;
    delete[] r_partial;
    delete[] l_sum;
    delete[] r_sum;
}

template<typename T>
void radix::OpenMPRadix<T>::operator ()(T* typed_a, size_t len) const
{
    typedef radix::radix_traits<T>::integer_type I;
    I *a = (I*)typed_a;
    I *d = new I[len];
    I* t;

    size_t work = len / num_threads + 1;;
    size_t i, mynum, last;

    for(I bit = 1; (bit & ~radix::radix_traits<T>::MSB_mask); bit <<=1, t = a, a = d, d = t)
    {
#pragma omp parallel default(none) private(i, mynum, last) shared(bit, a, d, t, work, len)
        {
            mynum = omp_get_thread_num();
            l_partial[mynum]= r_partial[mynum] = 0;

            for(i = work * mynum; i < work * mynum + work && i < len; i++)
               if(!(a[i] & bit))
                   l_partial[mynum] += 1;
               else
                   r_partial[mynum] += 1;
#pragma omp barrier

            l_sum[mynum]= r_sum[mynum] = 0;
            int border = 0;
            for(i = 0; i < num_threads; i++)
            {
                if(mynum > i)
                {
                    l_sum[mynum] += l_partial[i];
                    r_sum[mynum] += r_partial[i];
                }

                border += l_partial[i];
            }

#pragma omp barrier

            for(i = work * mynum; i < (last = work * mynum + work < len ? work * mynum + work : len); i++)
                if(!(a[i] & bit))
                    d[l_sum[mynum]++] = a[i];
                else
                    d[border + r_sum[mynum]++] = a[i];
        }
    }

    //MSD in inverse opder
#pragma omp parallel default(none) private(i, mynum, last) shared(a, d, t, work, len)
        {
            mynum = omp_get_thread_num();
            l_partial[mynum]= r_partial[mynum] = 0;

            for(i = work * mynum; i < work * mynum + work && i < len; i++)
               if((a[i] & radix::radix_traits<T>::MSB_mask))
                   l_partial[mynum] += 1;
               else
                   r_partial[mynum] += 1;
#pragma omp barrier
            l_sum[mynum]= r_sum[mynum] = 0;
            int border = 0;
            for(i = 0; i < num_threads; i++)
            {
                if(mynum > i)
                {
                    l_sum[mynum] += l_partial[i];
                    r_sum[mynum] += r_partial[i];
                }

                border += l_partial[i];
            }

#pragma omp barrier

            for(i = work * mynum; i < (last = work * mynum + work < len ? work * mynum + work : len); i++)
                if((a[i] & radix::radix_traits<T>::MSB_mask))
                    d[border -1 - l_sum[mynum]++] = a[i];
                else
                    d[border + r_sum[mynum]++] = a[i];
        }
    t = a; a = d; d = t;
    delete[] d;
}

void radix::OpenMPRadix<signed int>::operator ()(signed int* a, size_t len) const
{
    signed int *d = new signed int[len];
    signed int* t = 0;

    size_t work = len / num_threads + 1;
    size_t i, mynum, last;
    unsigned int bit = 1;
    for(bit = 1; (bit & ~0x80000000U) ; bit <<=1, t = a, a = d, d = t)
    {
#pragma omp parallel default(none) private(i, mynum, last) shared(bit, a, d, t, work, len)
        {
            mynum = omp_get_thread_num();
            l_partial[mynum]= r_partial[mynum] = 0;

            for(i = work * mynum; i < work * mynum + work && i < len; i++)
               if(!(a[i] & bit))
                   l_partial[mynum] += 1;
               else
                   r_partial[mynum] += 1;
#pragma omp barrier
            l_sum[mynum]= r_sum[mynum] = 0;
            int border = 0;
            for(i = 0; i < num_threads; i++)
            {
                if(mynum > i)
                {
                    l_sum[mynum] += l_partial[i];
                    r_sum[mynum] += r_partial[i];
                }

                border += l_partial[i];
            }

#pragma omp barrier

            for(i = work * mynum; i < (last = work * mynum + work < len ? work * mynum + work : len); i++)
                if(!(a[i] & bit))
                    d[l_sum[mynum]++] = a[i];
                else
                    d[border + r_sum[mynum]++] = a[i];
        }
    }
    //MSD

#pragma omp parallel default(none) private(i, mynum, last) shared(a, d, t, work, len)
        {
            mynum = omp_get_thread_num();
            l_partial[mynum]= r_partial[mynum] = 0;

            for(i = work * mynum; i < work * mynum + work && i < len; i++)
               if((a[i] & 0x80000000U))
                   l_partial[mynum] += 1;
               else
                   r_partial[mynum] += 1;
#pragma omp barrier
            l_sum[mynum]= r_sum[mynum] = 0;
            int border = 0;
            for(i = 0; i < num_threads; i++)
            {
                if(mynum > i)
                {
                    l_sum[mynum] += l_partial[i];
                    r_sum[mynum] += r_partial[i];
                }

                border += l_partial[i];
            }

#pragma omp barrier

            for(i = work * mynum; i < (last = work * mynum + work < len ? work * mynum + work : len); i++)
                if((a[i] & 0x80000000U))
                    d[l_sum[mynum]++] = a[i];
                else
                    d[border + r_sum[mynum]++] = a[i];
        }

    t = a; a = d; d = t;
    delete[] d;
}

template<>
void radix::OpenMPRadix<unsigned int>::operator ()(unsigned int* a, size_t len) const
{
    unsigned int *d = new unsigned int[len];
    unsigned int* t = 0;

    size_t work = len / num_threads + 1;
    size_t i, tid, last;
    unsigned int bit = 1;
    for(bit = 1; bit; bit <<=1, t = a, a = d, d = t)
    {
#pragma omp parallel default(none) private(i, tid, last) shared(bit, a, d, t, work, len)
        {
            tid = omp_get_thread_num();
            l_partial[tid]= r_partial[tid] = 0;

            for(i = work * tid; i < work * tid + work && i < len; i++)
               if(!(a[i] & bit))
                   l_partial[tid] += 1;
               else
                   r_partial[tid] += 1;
#pragma omp barrier
            l_sum[tid]= r_sum[tid] = 0;
            int border = 0;
            for(i = 0; i < num_threads; i++)
            {
                if(tid > i)
                {
                    l_sum[tid] += l_partial[i];
                    r_sum[tid] += r_partial[i];
                }

                border += l_partial[i];
            }

#pragma omp barrier

            for(i = work * tid; i < (last = work * tid + work < len ? work * tid + work : len); i++)
                if(!(a[i] & bit))
                    d[l_sum[tid]++] = a[i];
                else
                    d[border + r_sum[tid]++] = a[i];
        }
    }
    delete[] d;
}

template void radix::OpenMPRadix<float>::operator()(float*, size_t) const;
template void radix::OpenMPRadix<double>::operator()(double*, size_t) const;

template radix::OpenMPRadix<float>::OpenMPRadix();
template radix::OpenMPRadix<double>::OpenMPRadix();
template radix::OpenMPRadix<int>::OpenMPRadix();
template radix::OpenMPRadix<unsigned int>::OpenMPRadix();

template radix::OpenMPRadix<float>::~OpenMPRadix();
template radix::OpenMPRadix<double>::~OpenMPRadix();
template radix::OpenMPRadix<int>::~OpenMPRadix();
template radix::OpenMPRadix<unsigned int>::~OpenMPRadix();
