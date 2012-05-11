#include <sort_algorithm.hpp>
#include <traits.hpp>
#include <tbb/tbb.h>

template<typename T>
struct UnsignedOrder
{
    bool operator()(const T x, const T bit)
    {
        return !!(((x >> 1) ^ x) & bit);
    }
};

template<typename I, typename Ord>
class TBBRadixScanTask
{
public:
    TBBRadixScanTask(I* _a, I* _d, size_t _n, I _bit)
        : x(_a), y(_d), l_bound(0), r_bound(0), n(_n), bit(_bit){}

    template<typename Tag>
    void operator()( const tbb::blocked_range<int>& r, Tag )
    {
        Ord order;
        size_t old_l_bound = l_bound;
        size_t old_r_bound = r_bound;
        int i = 0;

        for(i = r.begin(); i < r.end(); ++i)
        {
            if (order(x[i], bit))
                ++r_bound;
            else
                ++l_bound;
        }

        if( Tag::is_final_scan() )
        {
            size_t k = 0;
            int m = 0;
            for(i = r.begin(); i < r.end(); ++i )
            {
                if (order(x[i], bit))
                    y[n - old_r_bound -1 + m--] = x[i];
                else
                    y[old_l_bound + k++] = x[i];
            }
        }

    }
    TBBRadixScanTask( TBBRadixScanTask& b, tbb::split ) : x(b.x), y(b.y),l_bound(0), r_bound(0), n(b.n), bit(b.bit) {}
    void reverse_join( TBBRadixScanTask& a )
    {
        l_bound = a.l_bound + l_bound;
        r_bound = a.r_bound + r_bound;
    }

    void assign( TBBRadixScanTask& b ) {l_bound = b.l_bound; r_bound = b.r_bound;}

private:
    const I* x;
    I* y;
    I bit;

    size_t l_bound;
    size_t r_bound;
    size_t n;
};

template<>
void radix::TBBRadix<unsigned int>::operator ()(unsigned int* a, size_t len) const
{
    unsigned int* d = new unsigned int[len];
    unsigned int* s = a;
    unsigned int* t;

    for (unsigned int bit = 1; bit; bit <<= 1, t = s, s = d, d = t)
    {
        TBBRadixScanTask<unsigned int, UnsignedOrder<unsigned int> > body(s, d, len, bit);
        tbb::parallel_scan( tbb::blocked_range<int>(0, len, 4), body);
    }
    delete[] d;
}

template<typename T>
void radix::TBBRadix<T>::operator ()(T* typed_a, size_t len) const
{
    typedef radix::radix_traits<T>::integer_type I;
    I *a = (I*)typed_a;
    I *d = new I[len];
    I* t;

    for (I bit = 1; (bit & ~radix::radix_traits<T>::MSB_mask); bit <<= 1, t = a, a = d, d = t)
    {
        TBBRadixScanTask<I, UnsignedOrder<I> > body(a, d, len, bit);
        tbb::parallel_scan( tbb::blocked_range<int>(0, len, 4), body);
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

template void radix::TBBRadix<float>::operator()(float*, size_t) const;
template void radix::TBBRadix<double>::operator()(double*, size_t) const;
template void radix::TBBRadix<signed int>::operator()(signed int*, size_t) const;
