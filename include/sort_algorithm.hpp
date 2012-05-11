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
    ~OpenMPRadix<T>();
    void operator ()(T* arr, size_t size) const;
private:
    size_t num_threads;
    unsigned int* partial;
    unsigned int* partial1;
    unsigned int* temp;
    unsigned int* temp1;
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
