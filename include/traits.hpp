#ifndef __traits_hpp__
#define __traits_hpp__

namespace radix
{
template<typename T>
struct radix_traits
{
    typedef T self_type;
};

template<>
struct radix_traits<signed int>
{
    typedef unsigned int integer_type;
    static const integer_type MSB_mask = 0x80000000U;

    static size_t index(size_t& l_bound, size_t middle)
    {
        return middle - l_bound;
    }
};

template<>
struct radix_traits<unsigned int>
{
    typedef unsigned int integer_type;
    static const integer_type MSB_mask = 0x80000000U;
};

template<>
struct radix_traits<float>
{
    typedef unsigned int integer_type;
    static const integer_type MSB_mask = 0x80000000U;

    static size_t index(size_t& l_bound, size_t middle)
    {
        return l_bound;
    }
};

template<>
struct radix_traits<double>
{
    typedef unsigned long long integer_type;
    static const integer_type MSB_mask = 0x8000000000000000ULL;

    static size_t index(size_t& l_bound, size_t middle)
    {
        return l_bound;
    }
};

}

#endif
