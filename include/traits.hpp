#ifndef __traits_hpp__
#define __traits_hpp__

#ifndef __CUDA_ARCH__
# include <vector>

namespace radix
{
class TestInfo
{
public:
    void add(std::pair<size_t, double> res)
    {
        measurements.push_back(res);
    }

public:
    std::vector<std::pair<size_t, double> > measurements;
};
}

#endif
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
    typedef unsigned int bit_extract_type;
    static const integer_type MSB_mask = 0x80000000U;
    static const integer_type MSB_bit = 31;

#ifdef __CUDA_ARCH__
    static __device__ __forceinline__ size_t index(size_t tid_offset, size_t middle)
    {
        if (!blockIdx.x & 1)
            return (tid_offset < middle) ? tid_offset + (blockDim.x << 1) - middle : tid_offset - middle;
        else
            return (tid_offset < middle) ? middle - 1 - tid_offset :((blockDim.x << 1) - (tid_offset - middle) - 1);
    }
    static __device__ __forceinline__ bit_extract_type as_integer(unsigned int x)
    {
        return x;
    }
#else
    static size_t index(size_t& l_bound, size_t middle)
    {
        return middle - l_bound;
    }
#endif
};

template<>
struct radix_traits<unsigned int>
{
    typedef unsigned int integer_type;
    typedef unsigned int bit_extract_type;
    static const integer_type MSB_mask = 0x80000000U;
    static const integer_type MSB_bit = 31;

#ifdef __CUDA_ARCH__
    static __device__ __forceinline__ size_t index(size_t tid_offset, size_t middle)
    {
        return (blockIdx.x & 1)? (blockDim.x << 1) - tid_offset - 1 : tid_offset;
    }
    static __device__ __forceinline__ bit_extract_type as_integer(signed int x)
    {
        return (bit_extract_type)x;
    }
#endif
};

template<>
struct radix_traits<float>
{
    typedef unsigned int integer_type;
    typedef signed int bit_extract_type;
    static const integer_type MSB_mask = 0x80000000U;
    static const integer_type MSB_bit = 31;

#ifdef __CUDA_ARCH__
    static __device__ __forceinline__ size_t index(size_t tid_offset, size_t middle)
    {
        if (!(blockIdx.x & 1))
            return (tid_offset < middle) ? tid_offset + (blockDim.x << 1) - middle : (blockDim.x << 1) - tid_offset -1;
        else
            return (tid_offset < middle) ? middle - 1 - tid_offset : tid_offset;
    }

    static __device__ __forceinline__ bit_extract_type as_integer(float x)
    {
        return __float_as_int(x);
    }
#else
    static size_t index(size_t& l_bound, size_t middle)
    {
        return l_bound;
    }
#endif
};

template<>
struct radix_traits<double>
{
    typedef unsigned long long integer_type;
    static const integer_type MSB_mask = 0x8000000000000000ULL;
    static const integer_type MSB_bit = 63;
    typedef signed long long bit_extract_type;
#ifdef __CUDA_ARCH__
    static __device__ __forceinline__ size_t index(size_t tid_offset, size_t middle)
    {
        if (!(blockIdx.x & 1))
            return (tid_offset < middle) ? tid_offset + (blockDim.x << 1) - middle : (blockDim.x << 1) - tid_offset -1;
        else
            return (tid_offset < middle) ? middle - 1 - tid_offset : tid_offset;
    }

    static __device__ __forceinline__ bit_extract_type as_integer(float x)
    {
        return __double_as_longlong(x);
    }
#else
    static size_t index(size_t& l_bound, size_t middle)
    {
        return l_bound;
    }
#endif
};

}

#endif
