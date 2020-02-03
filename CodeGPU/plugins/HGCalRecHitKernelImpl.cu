#include <cuda.h>
#include <cuda_runtime.h>
#include "HGCalRecHitKernelImpl.cuh"

__global__
void ee_step1(HGCUncalibratedRecHitSoA dst_soa, HGCUncalibratedRecHitSoA src_soa, const HGCeeUncalibratedRecHitConstantData cdata, size_t length)
{
  //shared memory allocation
  __shared__ double   sd[cdata.ndelem]; 
  __shared__ float    sf[cdata.nfelem]; 
  __shared__ uint32_t su[cdata.nuelem];
  __shared__ bool     sb[cdata.nbelem];
  
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

  //setting shared memory
  sd[0] = cdata.hgcEE_keV2DIGI_;
  sd[1] = cdata.hgceeUncalib2GeV_;
  for(unsigned int i=0 && i < cdata.s_hgcEE_fCPerMIP_ )
    sd[i+2] = cdata.hgcEE_fCPerMIP_[i];
  for(unsigned int i=0 && i < cdata.s_hgcEE_cce_ )
    sd[i+2+cdata.s_hgcEE_fCPerMIP_] = cdata.hgcEE_cce_[i];
  for(unsigned int i=0 && i < cdata.s_hgcEE_noise_fC_ )
    sd[i+2+cdata.s_hgcEE_fCPerMIP_+cdata.s_hgcEE_cce_] = cdata.hgcEE_noise_fC_[i];
  for(unsigned int i=0 && i < cdata.s_rcorr_ )
    sd[i+2+cdata.s_hgcEE_fCPerMIP_+cdata.s_hgcEE_cce_+cdata.s_hgcEE_noise_fC_] = cdata.rcorr_[i];

  for(unsigned int i=0 && i < cdata.s_weights_ )
    sf[i] = cdata.weights_[i];

  su[0] = cdata.rangeMatch_;
  su[1] = cdata.rangeMask_;

  sb[0] = hgcEE_isSiFE_;

  __syncthreads();

  for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < length; i += blockDim.x * gridDim.x)
    {
      dst_soa.amplitude[i] = src_soa.amplitude[i];
    }
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
