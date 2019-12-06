#ifndef _HGCALRECHITKERNELIMPL_CUH_
#define _HGCALRECHITKERNELIMPL_CUH_

#include <cuda.h>
#include <cuda_runtime.h>

class HGCRecHit_GPU {
  float energy;
  float time;
  typedef int key_type;
};

class HGCUncalibratedRecHit_GPU {
  float amplitude;
  float pedestal;
  typedef int key_type;
};

/*
__global__
void ee_step1(HGCUncalibratedRecHit *__restrict__ dst, HGCUncalibratedRecHit *__restrict__ src, size_t length);
__global__
void hef_step1(HGCUncalibratedRecHit *__restrict__ dst, HGCUncalibratedRecHit *__restrict__ src, size_t length);
__global__
void heb_step1(HGCUncalibratedRecHit *__restrict__ dst, HGCUncalibratedRecHit *__restrict__ src, size_t length);
__global__
void to_rechit(HGCRecHit *__restrict__ dst, HGCUncalibratedRecHit *__restrict__ src, size_t length);
*/
__global__
void ee_step1(HGCUncalibratedRecHit_GPU *__restrict__ dst, HGCUncalibratedRecHit_GPU *__restrict__ src, size_t length);
__global__
void hef_step1(HGCUncalibratedRecHit_GPU *__restrict__ dst, HGCUncalibratedRecHit_GPU *__restrict__ src, size_t length);
__global__
void heb_step1(HGCUncalibratedRecHit_GPU *__restrict__ dst, HGCUncalibratedRecHit_GPU *__restrict__ src, size_t length);
__global__
void to_rechit(HGCRecHit_GPU *__restrict__ dst, HGCUncalibratedRecHit_GPU *__restrict__ src, size_t length);

#endif /* _HGCALRECHITKERNELIMPL_H_ */
