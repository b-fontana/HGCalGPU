#ifndef _HGCALRECHITKERNELIMPL_CUH_
#define _HGCALRECHITKERNELIMPL_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include "Types.h"

/*
=======
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
void ee_step1(HGCUncalibratedRecHitSoA *__restrict__ dst, HGCUncalibratedRecHitSoA *__restrict__ src, size_t length);
__global__
void hef_step1(HGCUncalibratedRecHitSoA *__restrict__ dst, HGCUncalibratedRecHitSoA *__restrict__ src, size_t length);
__global__
void heb_step1(HGCUncalibratedRecHitSoA *__restrict__ dst, HGCUncalibratedRecHitSoA *__restrict__ src, size_t length);
__global__
void to_rechit(HGCRecHitSoA *__restrict__ dst, HGCUncalibratedRecHitSoA *__restrict__ src, size_t length);

#endif /* _HGCALRECHITKERNELIMPL_H_ */
