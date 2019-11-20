#ifndef HeterogeneousUncalibRecHitsProd_cuh
#define HeterogeneousUncalibRecHitsProd_cuh

#include <cuda_runtime.h>

__global__ 
void scrambler_kernel(const char * message, size_t length);

void scrambler_wrapper(const char * message, size_t length);

#endif //HeterogeneousUncalibRecHitsProd_cuh
