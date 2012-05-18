#ifndef __cuda_radix_h__
#define __cuda_radix_h__

void cuda_radix_coller(unsigned int* a, unsigned int* d, size_t len);
void cuda_radix_coller(signed int* a, signed int* d, size_t len);
void cuda_radix_coller(float* a, float* d, size_t len);
void cuda_radix_coller(double* a, double* d, size_t len);

#endif
