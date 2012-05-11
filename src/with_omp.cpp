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
    partial =   new unsigned int[num_threads];
    partial1 =  new unsigned int[num_threads];
    temp =      new unsigned int[num_threads];
    temp1 =     new unsigned int[num_threads];
}

template<typename T>
radix::OpenMPRadix<T>::~OpenMPRadix()
{
    delete[] partial;
    delete[] partial1;
    delete[] temp;
    delete[] temp1;
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
            partial[mynum]= partial1[mynum] = 0;

            for(i = work * mynum; i < work * mynum + work && i < len; i++)
               if(!(a[i] & bit))
                   partial[mynum] += 1;
               else
                   partial1[mynum] += 1;
#pragma omp barrier

            temp[mynum]= temp1[mynum] = 0;
            int border = 0;
            for(i = 0; i < num_threads; i++)
            {
                if(mynum > i)
                {
                    temp[mynum] += partial[i];
                    temp1[mynum] += partial1[i];
                }

                border += partial[i];
            }

#pragma omp barrier

            for(i = work * mynum; i < (last = work * mynum + work < len ? work * mynum + work : len); i++)
                if(!(a[i] & bit))
                    d[temp[mynum]++] = a[i];
                else
                    d[border + temp1[mynum]++] = a[i];
        }
    }

    //MSD in inverse opder
#pragma omp parallel default(none) private(i, mynum, last) shared(a, d, t, work, len)
        {
            mynum = omp_get_thread_num();
            partial[mynum]= partial1[mynum] = 0;

            for(i = work * mynum; i < work * mynum + work && i < len; i++)
               if((a[i] & radix::radix_traits<T>::MSB_mask))
                   partial[mynum] += 1;
               else
                   partial1[mynum] += 1;
#pragma omp barrier
            temp[mynum]= temp1[mynum] = 0;
            int border = 0;
            for(i = 0; i < num_threads; i++)
            {
                if(mynum > i)
                {
                    temp[mynum] += partial[i];
                    temp1[mynum] += partial1[i];
                }

                border += partial[i];
            }

#pragma omp barrier

            for(i = work * mynum; i < (last = work * mynum + work < len ? work * mynum + work : len); i++)
                if((a[i] & radix::radix_traits<T>::MSB_mask))
                    d[border -1 - temp[mynum]++] = a[i];
                else
                    d[border + temp1[mynum]++] = a[i];
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
            partial[mynum]= partial1[mynum] = 0;

            for(i = work * mynum; i < work * mynum + work && i < len; i++)
               if(!(a[i] & bit))
                   partial[mynum] += 1;
               else
                   partial1[mynum] += 1;
#pragma omp barrier
            temp[mynum]= temp1[mynum] = 0;
            int border = 0;
            for(i = 0; i < num_threads; i++)
            {
                if(mynum > i)
                {
                    temp[mynum] += partial[i];
                    temp1[mynum] += partial1[i];
                }

                border += partial[i];
            }

#pragma omp barrier

            for(i = work * mynum; i < (last = work * mynum + work < len ? work * mynum + work : len); i++)
                if(!(a[i] & bit))
                    d[temp[mynum]++] = a[i];
                else
                    d[border + temp1[mynum]++] = a[i];
        }
    }
    //MSD

#pragma omp parallel default(none) private(i, mynum, last) shared(a, d, t, work, len)
        {
            mynum = omp_get_thread_num();
            partial[mynum]= partial1[mynum] = 0;

            for(i = work * mynum; i < work * mynum + work && i < len; i++)
               if((a[i] & 0x80000000U))
                   partial[mynum] += 1;
               else
                   partial1[mynum] += 1;
#pragma omp barrier
            temp[mynum]= temp1[mynum] = 0;
            int border = 0;
            for(i = 0; i < num_threads; i++)
            {
                if(mynum > i)
                {
                    temp[mynum] += partial[i];
                    temp1[mynum] += partial1[i];
                }

                border += partial[i];
            }

#pragma omp barrier

            for(i = work * mynum; i < (last = work * mynum + work < len ? work * mynum + work : len); i++)
                if((a[i] & 0x80000000U))
                    d[temp[mynum]++] = a[i];
                else
                    d[border + temp1[mynum]++] = a[i];
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
    size_t i, mynum, last;
    unsigned int bit = 1;
    for(bit = 1; bit; bit <<=1, t = a, a = d, d = t)
    {
#pragma omp parallel default(none) private(i, mynum, last) shared(bit, a, d, t, work, len)
        {
            mynum = omp_get_thread_num();
            partial[mynum]= partial1[mynum] = 0;

            for(i = work * mynum; i < work * mynum + work && i < len; i++)
               if(!(a[i] & bit))
                   partial[mynum] += 1;
               else
                   partial1[mynum] += 1;
#pragma omp barrier
            temp[mynum]= temp1[mynum] = 0;
            int border = 0;
            for(i = 0; i < num_threads; i++)
            {
                if(mynum > i)
                {
                    temp[mynum] += partial[i];
                    temp1[mynum] += partial1[i];
                }

                border += partial[i];
            }

#pragma omp barrier

            for(i = work * mynum; i < (last = work * mynum + work < len ? work * mynum + work : len); i++)
                if(!(a[i] & bit))
                    d[temp[mynum]++] = a[i];
                else
                    d[border + temp1[mynum]++] = a[i];
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
