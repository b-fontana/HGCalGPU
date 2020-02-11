#include <cuda.h>
#include <cuda_runtime.h>
#include "HGCalRecHitKernelImpl.cuh"

__device__
int wafer(uint32_t id)
{
  static const int kHGCalWaferOffset = 8;
  static const int kHGCalWaferMask = 0x3FF;
  return (id >> kHGCalWaferOffset) & kHGCalWaferMask; 
}

__device__ 
void set_shared_memory(int tid, double*& sd, float*& sf, uint32_t*& su, int*& si, bool*& sb, const HGCeeUncalibratedRecHitConstantData& cdata)
{
  size_t size1 = cdata.s_hgcEE_fCPerMIP_ + 2;
  size_t size2 = cdata.s_hgcEE_cce_      + size1;
  size_t size3 = cdata.s_hgcEE_noise_fC_ + size2;
  size_t size4 = cdata.s_rcorr_          + size3; 
  size_t size5 = cdata.s_weights_        + size4; 
  size_t size6 = cdata.s_waferTypeL_     + size5; 
  if(tid == 0)
    printf("%f, %f, %f, %f, %f, %f", size1, size2, size3, size4, size5, size6);

  if(tid == 0)
    sd[tid] = cdata.hgcEE_keV2DIGI_;
  else if(tid == 1)
    sd[tid] = cdata.hgceeUncalib2GeV_;
  else if(tid > 1 && tid < size1)
    sd[tid] = cdata.hgcEE_fCPerMIP_[tid-2];
  else if(tid >= size1 && tid < size2)
    sd[tid] = cdata.hgcEE_cce_[tid-size1];
  else if(tid >= size2 && tid < size3)
    sd[tid] = cdata.hgcEE_noise_fC_[tid-size2];
  else if(tid >= size3 && tid < size4)
    sd[tid] = cdata.rcorr_[tid - size3];
  else if(tid >= size4 && tid < size5)
    sd[tid] = cdata.weights_[tid - size4];
  else if(tid >= size5 && tid < size6)
    si[tid] = cdata.waferTypeL_[tid - size5];
  else if(tid == size6)
    su[tid] = cdata.rangeMatch_;
  else if(tid == size6 + 1)
    su[tid] = cdata.rangeMask_;
  else if(tid == size6 + 2)
    sb[tid] = cdata.hgcEE_isSiFE_;

  __syncthreads();
}

__global__
void ee_step1(HGCUncalibratedRecHitSoA dst_soa, HGCUncalibratedRecHitSoA src_soa, const HGCeeUncalibratedRecHitConstantData cdata, size_t length)
{
  extern __shared__ double s[];
  /*
  printf("check\n");
  __syncthreads();
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if(tid==0)
    {
      printf("%zu, %zu, %zu, %zu, %zu", cdata.ndelem, cdata.nfelem, cdata.nuelem, cdata.nielem, cdata.nbelem);

      extern __shared__ double s[];
      double   *sd = s;
      float    *sf = (float*)   (sd + cdata.ndelem);
      uint32_t *su = (uint32_t*)(sf + cdata.nfelem);
      int      *si = (int*)     (su + cdata.nuelem);
      printf("here 5\n");
      bool     *sb = (bool*)    (si + cdata.nielem); //erro!!!!!
      printf("here 6\n");
      set_shared_memory(threadIdx.x, sd, sf, su, si, sb, cdata);
      printf("here 7\n");

      printf("%f\n", sd[2]);
      printf("%d\n", si[3]);
      printf("%d\n", su[1]);
      printf("%d\n", sb[0]);
    }

  __syncthreads();

  for (unsigned int i = tid; i < length; i += blockDim.x * gridDim.x)
    {
      dst_soa.amplitude[i] = src_soa.amplitude[i];
    }
  */
}

__global__
void hef_step1(HGCUncalibratedRecHitSoA dst_soa, HGCUncalibratedRecHitSoA src_soa, const HGChefUncalibratedRecHitConstantData cdata, size_t length)
{
}

__global__
void heb_step1(HGCUncalibratedRecHitSoA dst_soa, HGCUncalibratedRecHitSoA src_soa, const HGChebUncalibratedRecHitConstantData cdata, size_t length)
{
}

__global__
void ee_to_rechit(HGCRecHitSoA dst_soa, HGCUncalibratedRecHitSoA src_soa, const HGCeeUncalibratedRecHitConstantData cdata, size_t length)
{
  for (size_t i = blockDim.x * blockIdx.x + threadIdx.x; i < length; i += blockDim.x * gridDim.x)
    {
      dst_soa.energy[i] = 2.;
    }
}

__global__
void hef_to_rechit(HGCRecHitSoA dst_soa, HGCUncalibratedRecHitSoA src_soa, const HGChefUncalibratedRecHitConstantData cdata, size_t length)
{
}

__global__
void heb_to_rechit(HGCRecHitSoA dst_soa, HGCUncalibratedRecHitSoA src_soa, const HGChebUncalibratedRecHitConstantData cdata, size_t length)
{
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
