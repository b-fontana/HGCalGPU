#ifndef _HGCALRECHITKERNELIMPL_CUH_
#define _HGCALRECHITKERNELIMPL_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include "Types.h"

__global__
void ee_step1(float *__restrict__ dst_amplitude, float *__restrict__ src_amplitude, float *__restrict__ dst_pedestal, float *__restrict__ src_pedestal, float *__restrict__ dst_jitter, float *__restrict__ src_jitter, float *__restrict__ dst_chi2, float *__restrict__ src_chi2, float *__restrict__ dst_OOTamplitude, float *__restrict__ src_OOTamplitude, float *__restrict__ dst_OOTchi2, float *__restrict__ src_OOTchi2, uint32_t *__restrict__ dst_flags, uint32_t *__restrict__ src_flags, uint32_t *__restrict__ dst_aux, uint32_t *__restrict__ src_aux, uint32_t *__restrict__ dst_id, uint32_t *__restrict__ src_id, size_t length);
__global__
void ee_step2(float *__restrict__ dst_amplitude, float *__restrict__ src_amplitude, float *__restrict__ dst_pedestal, float *__restrict__ src_pedestal, float *__restrict__ dst_jitter, float *__restrict__ src_jitter, float *__restrict__ dst_chi2, float *__restrict__ src_chi2, float *__restrict__ dst_OOTamplitude, float *__restrict__ src_OOTamplitude, float *__restrict__ dst_OOTchi2, float *__restrict__ src_OOTchi2, uint32_t *__restrict__ dst_flags, uint32_t *__restrict__ src_flags, uint32_t *__restrict__ dst_aux, uint32_t *__restrict__ src_aux, uint32_t *__restrict__ dst_id, uint32_t *__restrict__ src_id, size_t length);
__global__
void hef_step1(HGCUncalibratedRecHitSoA *__restrict__ dst, HGCUncalibratedRecHitSoA *__restrict__ src, size_t length);
__global__
void heb_step1(HGCUncalibratedRecHitSoA *__restrict__ dst, HGCUncalibratedRecHitSoA *__restrict__ src, size_t length);
__global__
void to_rechit(HGCRecHitSoA *__restrict__ dst, HGCUncalibratedRecHitSoA *__restrict__ src, size_t length);

#endif /* _HGCALRECHITKERNELIMPL_H_ */
