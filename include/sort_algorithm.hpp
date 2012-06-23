#ifndef __sort_algorithm_h__
#define __sort_algorithm_h__
#include <traits.hpp>
//#include <algorithm>

namespace radix
{
template<typename T>
class SequentalSTD
{
public:
    SequentalSTD(){}
    void operator ()(T* arr, size_t size)
    {
        std::sort(arr, arr+ size);
    }
};

template<typename T>
class SequentalRadix
{
public:
    SequentalRadix(){}
    double operator ()(T* arr, size_t size) const;
};

template<typename T>
class TBBRadix
{
public:
    TBBRadix(){}
    double operator ()(T* arr, size_t size) const;
};

template<typename T>
class OpenMPRadix
{
public:
    OpenMPRadix();
    ~OpenMPRadix();
    double operator ()(T* arr, size_t size) const;
private:
    size_t num_threads;
    unsigned int* l_partial;
    unsigned int* r_partial;
    unsigned int* l_sum;
    unsigned int* r_sum;
};

template<typename T>
class CudaRadix
{
public:
    CudaRadix(){}
    double operator ()(T* arr, size_t size) const;
};

template<typename T>
class SilkRadix
{
public:
    SilkRadix(){}
    double operator ()(T* arr, size_t size) const;
};
}

#endif
