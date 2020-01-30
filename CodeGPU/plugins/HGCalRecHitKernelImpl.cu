#include <cuda.h>
#include <cuda_runtime.h>
#include "HGCalRecHitKernelImpl.cuh"

__global__
void ee_step1(float *__restrict__ dst_amplitude, float *__restrict__ dst_pedestal, float *__restrict__ dst_jitter,  float *__restrict__ dst_chi2, float *__restrict__ dst_OOTamplitude, float *__restrict__ dst_OOTchi2, uint32_t *__restrict__ dst_flags, uint32_t *__restrict__ dst_aux, uint32_t *__restrict__ dst_id, float *__restrict__ src_amplitude, float *__restrict__ src_pedestal, float *__restrict__ src_jitter,  float *__restrict__ src_chi2, float *__restrict__ src_OOTamplitude, float *__restrict__ src_OOTchi2, uint32_t *__restrict__ src_flags, uint32_t *__restrict__ src_aux, uint32_t *__restrict__ src_id, double cdata, size_t length)
{
  if(blockDim.x * blockIdx.x + threadIdx.x == 0)
    printf("Constant data: %f\n", cdata);
  for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < length; i += blockDim.x * gridDim.x)
    {
      dst_amplitude[i] = src_amplitude[i];
    }
}

__global__
void hef_step1(float *__restrict__ dst_amplitude, float *__restrict__ dst_pedestal, float *__restrict__ dst_jitter,  float *__restrict__ dst_chi2, float *__restrict__ dst_OOTamplitude, float *__restrict__ dst_OOTchi2, uint32_t *__restrict__ dst_flags, uint32_t *__restrict__ dst_aux, uint32_t *__restrict__ dst_id, float *__restrict__ src_amplitude, float *__restrict__ src_pedestal, float *__restrict__ src_jitter,  float *__restrict__ src_chi2, float *__restrict__ src_OOTamplitude, float *__restrict__ src_OOTchi2, uint32_t *__restrict__ src_flags, uint32_t *__restrict__ src_aux, uint32_t *__restrict__ src_id, size_t length)
{
}

__global__
void heb_step1(float *__restrict__ dst_amplitude, float *__restrict__ dst_pedestal, float *__restrict__ dst_jitter,  float *__restrict__ dst_chi2, float *__restrict__ dst_OOTamplitude, float *__restrict__ dst_OOTchi2, uint32_t *__restrict__ dst_flags, uint32_t *__restrict__ dst_aux, uint32_t *__restrict__ dst_id, float *__restrict__ src_amplitude, float *__restrict__ src_pedestal, float *__restrict__ src_jitter,  float *__restrict__ src_chi2, float *__restrict__ src_OOTamplitude, float *__restrict__ src_OOTchi2, uint32_t *__restrict__ src_flags, uint32_t *__restrict__ src_aux, uint32_t *__restrict__ src_id, size_t length)
{
}

__global__
void to_rechit(float *__restrict__ dst_energy, float *__restrict__ dst_time, uint32_t *__restrict__ dst_id, uint32_t *__restrict__ dst_flags, uint32_t *__restrict__ dst_flagBits, float *__restrict__ src_amplitude, float *__restrict__ src_pedestal, float *__restrict__ src_jitter,  float *__restrict__ src_chi2, float *__restrict__ src_OOTamplitude, float *__restrict__ src_OOTchi2, uint32_t *__restrict__ src_flags, uint32_t *__restrict__ src_aux, uint32_t *__restrict__ src_id, size_t length)
{
  for (size_t i = blockDim.x * blockIdx.x + threadIdx.x; i < length; i += blockDim.x * gridDim.x)
    {
      dst_energy[i] = 2.;
    }
}

/*
=======
>>>>>>> b5bfc7e2f47f926abb3dcd21cdf5e2094e53dd3f
//declared as extern in DataFormats/CaloRecHit/interface/CaloRecHit.h
#ifdef __CUDA_ARCH__
__constant__ uint32_t calo_rechit_masks[] = {0x00000000u, 0x00000001u, 0x00000003u, 0x00000007u, 0x0000000fu, 0x0000001fu,
					     0x0000003fu, 0x0000007fu, 0x000000ffu, 0x000001ffu, 0x000003ffu, 0x000007ffu,
					     0x00000fffu, 0x00001fffu, 0x00003fffu, 0x00007fffu, 0x0000ffffu, 0x0001ffffu,
					     0x0003ffffu, 0x0007ffffu, 0x000fffffu, 0x001fffffu, 0x003fffffu, 0x007fffffu,
					     0x00ffffffu, 0x01ffffffu, 0x03ffffffu, 0x07ffffffu, 0x0fffffffu, 0x1fffffffu,
					     0x3fffffffu, 0x7fffffffu, 0xffffffffu};
#endif
<<<<<<< HEAD
*/
