#ifndef _KERNELMANAGER_H
#define _KERNELMANAGER_H

#include "DataFormats/DetId/interface/DetId.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCompat.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HGCalRecHitKernelImpl.cuh"
#include "Types.h"

#include <vector>
#include <algorithm> //std::swap  
#include <variant>
#include <cuda.h>
#include <cuda_runtime.h>

/*
#ifdef __CUDA_ARCH__
extern __constant__ uint32_t calo_rechit_masks[];
#endif
*/

template <typename TYPE_IN, typename TYPE_OUT>
  class KernelConstantData {
 public:
 KernelConstantData(double hgcEE_keV2DIGI_): hgcEE_keV2DIGI(hgcEE_keV2DIGI_) {}
  
  double hgcEE_keV2DIGI;
};

template <typename TYPE_IN, typename TYPE_OUT>
  class KernelModifiableData {
 public:
 KernelModifiableData(size_t nhits_, TYPE_IN *h_in_, TYPE_IN *d_1_, TYPE_IN *d_2_, TYPE_OUT *d_out_, TYPE_OUT *h_out_):
  nhits(nhits_), h_in(h_in_), d_1(d_1_), d_2(d_2_), d_out(d_out_), h_out(h_out_) {}

  size_t nhits;
  TYPE_IN *h_in;
  TYPE_IN *d_1, *d_2;
  TYPE_OUT *d_out;
  TYPE_OUT *h_out;
};

class KernelManagerBase {
public:
  explicit KernelManagerBase(){};
  virtual ~KernelManagerBase(){};
  
  virtual void run_kernels() = 0;

protected:
  virtual void assign_and_transfer_to_device() = 0;
  virtual void transfer_to_host_and_synchronize() = 0;
  virtual void reuse_device_pointers() = 0;
};

//the class assumes that the sizes of the arrays pointed to and the size of the collection are constant
class KernelManagerHGCalRecHit: private KernelManagerBase {
 public:
  explicit KernelManagerHGCalRecHit(KernelModifiableData<HGCUncalibratedRecHitSoA, HGCRecHitSoA>, const KernelConstantData<HGCUncalibratedRecHitSoA, HGCRecHitSoA>&, const DetId::Detector&);
  ~KernelManagerHGCalRecHit();
  void run_kernels();
  HGCRecHitSoA* get_output();

 private:
  void ee_step1_wrapper();
  void hef_step1_wrapper();
  void heb_step1_wrapper();
  void to_rechit_wrapper();
  void assign_and_transfer_to_device() override;
  void transfer_to_host_and_synchronize() override;
  void reuse_device_pointers() override;

  const DetId::Detector dtype_;
  const std::vector<HGCUncalibratedRecHitSoA> h_oldhits_;
  KernelModifiableData<HGCUncalibratedRecHitSoA, HGCRecHitSoA> data_;
  KernelConstantData<HGCUncalibratedRecHitSoA, HGCRecHitSoA> cdata_;
};

#endif //_KERNELMANAGER_H_
