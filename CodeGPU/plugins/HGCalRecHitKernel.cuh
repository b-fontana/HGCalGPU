#ifndef HGCalRecHitKernel_cuh
#define HGCalRecHitKernel_cuh

#include <cuda_runtime.h>

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

__global__ void scrambler_kernel(const char * message, size_t length);

__host__ void scrambler_wrapper(const char * message, size_t length);

__global__ void HGCeeRecHitKernel(HGCUncalibratedRecHit*, HGCUncalibratedRecHit*, size_t);

__host__ std::unique_ptr<HGCeeUncalibratedRecHitCollection> HGCeeRecHitKernel_wrapper(HGCeeUncalibratedRecHitCollection, 
										      const std::string&);

#endif //HGCalRecHitKernel_cuh
