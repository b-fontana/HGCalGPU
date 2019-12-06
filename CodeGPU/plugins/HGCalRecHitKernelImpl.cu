#include <cuda.h>
#include <cuda_runtime.h>
#include <inttypes.h>
#include "HGCalRecHitKernelImpl.cuh"

__device__
int wafer(uint32_t id)
{
  static const int kHGCalWaferOffset = 8;
  static const int kHGCalWaferMask = 0x3FF;
  return (id >> kHGCalWaferOffset) & kHGCalWaferMask; 
}

__device__
int layer(uint32_t id)
{
  static const int kHGCalLayerOffset = 20;
  static const int kHGCalLayerMask = 0x1F;
  return (id >> kHGCalLayerOffset) & kHGCalLayerMask; 
}

__device__ 
double get_weight_from_layer(const LENGTHSIZE& padding, const int& layer, double*& sd)
{
  return sd[padding + layer];
}

__device__ 
double get_thickness_correction(const LENGTHSIZE& padding, double *& sd, const HGCeeUncalibratedRecHitConstantData& cdata)
{
  int waferTypeL = cdata.waferTypeL_[1]; //this should be obtained fro the DetId and the wafer() __device__ function.
  return sd[padding + waferTypeL];
}

__device__
double get_noise(const LENGTHSIZE& padding, double *& sd, const HGCeeUncalibratedRecHitConstantData& cdata)
{
  int waferTypeL = cdata.waferTypeL_[1]; //this should be obtained fro the DetId and the wafer() __device__ function.
  return sd[padding + waferTypeL - 1];
}

__device__
double get_cce_correction(const LENGTHSIZE& padding, double *& sd, const HGCeeUncalibratedRecHitConstantData& cdata)
{
  int waferTypeL = cdata.waferTypeL_[1]; //this should be obtained fro the DetId and the wafer() __device__ function.
  return sd[padding + waferTypeL - 1];
}

__device__ 
double get_fCPerMIP(const LENGTHSIZE& padding, double *& sd, const HGCeeUncalibratedRecHitConstantData& cdata)
{
  int waferTypeL = cdata.waferTypeL_[1]; //this should be obtained fro the DetId and the wafer() __device__ function.
  return sd[padding + waferTypeL - 1];
}

__device__ 
void set_shared_memory(int tid, double*& sd, float*& sf, uint32_t*& su, int*& si, bool*& sb, const HGCeeUncalibratedRecHitConstantData& cdata, const LENGTHSIZE& size1, const LENGTHSIZE& size2, const LENGTHSIZE& size3, const LENGTHSIZE& size4, const LENGTHSIZE& size5, const LENGTHSIZE& size6)
{
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
    si[tid - size5] = cdata.waferTypeL_[tid - size5];
  else if(tid == size6)
    su[0] = cdata.rangeMatch_;
  else if(tid == size6 + 1)
    su[1] = cdata.rangeMask_;
  else if(tid == size6 + 2)
    sb[0] = cdata.hgcEE_isSiFE_;

  __syncthreads();
}

__global__
void ee_step1(HGCUncalibratedRecHitSoA dst_soa, HGCUncalibratedRecHitSoA src_soa, const HGCeeUncalibratedRecHitConstantData cdata, LENGTHSIZE length)
{
  //this kernel is currently doing nothing
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
  for (unsigned int i = tid; i < length; i += blockDim.x * gridDim.x)
    {
      printf("QQQQ %" PRIu32 "\n", src_soa.id[i]);
    }
}

__global__
void hef_step1(HGCUncalibratedRecHitSoA dst_soa, HGCUncalibratedRecHitSoA src_soa, const HGChefUncalibratedRecHitConstantData cdata, LENGTHSIZE length)
{
}

__global__
void heb_step1(HGCUncalibratedRecHitSoA dst_soa, HGCUncalibratedRecHitSoA src_soa, const HGChebUncalibratedRecHitConstantData cdata, LENGTHSIZE length)
{
}

__global__
void ee_to_rechit(HGCRecHitSoA dst_soa, HGCUncalibratedRecHitSoA src_soa, const HGCeeUncalibratedRecHitConstantData cdata, LENGTHSIZE length)
{
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

  int size1 = cdata.s_hgcEE_fCPerMIP_ + 2;
  int size2 = cdata.s_hgcEE_cce_      + size1;
  int size3 = cdata.s_hgcEE_noise_fC_ + size2;
  int size4 = cdata.s_rcorr_          + size3; 
  int size5 = cdata.s_weights_        + size4; 
  int size6 = cdata.s_waferTypeL_     + size5; 

  extern __shared__ double s[];
  double   *sd = s;
  float    *sf = (float*)   (sd + cdata.ndelem);
  uint32_t *su = (uint32_t*)(sf + cdata.nfelem);
  int      *si = (int*)     (su + cdata.nuelem);
  bool     *sb = (bool*)    (si + cdata.nielem);
  set_shared_memory(threadIdx.x, sd, sf, su, si, sb, cdata, size1, size2, size3, size4, size5, size6);

  for (unsigned int i = tid; i < length; i += blockDim.x * gridDim.x)
    {
      dst_soa.id[i] = src_soa.id[i];
      printf("ID: %" PRIu32 ", %" PRIu32 "\n", src_soa.id[i], dst_soa.id[i]);
      int l = layer(dst_soa.id[i]);
      double weight = get_weight_from_layer(size4, l, sd);
      double rcorr = get_thickness_correction(size3, sd, cdata);
      double noise = get_noise(size2, sd, cdata);
      double cce_correction = get_cce_correction(size1, sd, cdata);
      double fCPerMIP = get_fCPerMIP(2, sd, cdata);
      double sigmaNoiseGeV = 1e-3 * weight * rcorr * noise / fCPerMIP;

      //makeRecHit
      dst_soa.energy[i] = src_soa.amplitude[i] * weight * 0.001f * rcorr / cce_correction;
      dst_soa.time[i] = src_soa.jitter[i];
      dst_soa.flagBits[i] |= (0x1 << FlagsGPU::kGood);
      dst_soa.son[i] = nearbyintf( fmaxf(32.f, dst_soa.energy[i]/sigmaNoiseGeV) / ( 32.f * ((1 << 8)-1) ) );
      dst_soa.timeError[i] = 1.f; //change!!!!!!!!!!!!

      printf("quant %d: %d, %d, %f, %f, %f, %f, %f, %f\n", i, l, dst_soa.id[i], weight, rcorr, noise, cce_correction, fCPerMIP, sigmaNoiseGeV);
      printf("src %d: %f, %f, %f, %f, %f, %f, %d, %d, %d\n", i, src_soa.amplitude[i], src_soa.pedestal[i], src_soa.jitter[i], src_soa.chi2[i], src_soa.OOTamplitude[i], src_soa.OOTchi2[i], (int)src_soa.flags[i], (int)src_soa.aux[i], (int)src_soa.id[i]);
      printf("dst %d: %f, %f, %f, %d, %d, %d\n", i, dst_soa.energy[i], dst_soa.time[i], dst_soa.timeError[i], dst_soa.id[i], dst_soa.flagBits[i], dst_soa.son[i]);
    }
}

__global__
void hef_to_rechit(HGCRecHitSoA dst_soa, HGCUncalibratedRecHitSoA src_soa, const HGChefUncalibratedRecHitConstantData cdata, LENGTHSIZE length)
{
}

__global__
void heb_to_rechit(HGCRecHitSoA dst_soa, HGCUncalibratedRecHitSoA src_soa, const HGChebUncalibratedRecHitConstantData cdata, LENGTHSIZE length)
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
