#ifndef __sort_algorithm_h__
#define __sort_algorithm_h__
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
    void operator ()(T* arr, size_t size) const;
};

template<typename T>
class TBBRadix
{
public:
    TBBRadix(){}
    void operator ()(T* arr, size_t size) const;
};

template<typename T>
class OpenMPRadix
{
public:
    OpenMPRadix();
    ~OpenMPRadix();
    void operator ()(T* arr, size_t size) const;
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
    void operator ()(T* arr, size_t size) const;
};

template<typename T>
class SilkRadix
{
public:
    SilkRadix(){}
    void operator ()(T* arr, size_t size) const;
};
}

#endif
