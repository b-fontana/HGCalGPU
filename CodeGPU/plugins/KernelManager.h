#ifndef _KERNELMANAGER_H
#define _KERNELMANAGER_H

#include "FWCore/Utilities/interface/Exception.h"
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

template <typename T>
class KernelConstantData {
 public:
 KernelConstantData(const T& data_): data(data_) {
    if( std::is_same<T, HGCeeUncalibratedRecHitConstantData>::value or	std::is_same<T, HGChefUncalibratedRecHitConstantData>::value or std::is_same<T, HGChebUncalibratedRecHitConstantData>::value )
      {
	throw cms::Exception("WrongTemplateType") << "The KernelConstantData class does not support this type.";
      }
  }
  T data;
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

class KernelManagerHGCalRecHit {
 public:
  KernelManagerHGCalRecHit(KernelModifiableData<HGCUncalibratedRecHitSoA, HGCRecHitSoA>, const DetId::Detector&);
  ~KernelManagerHGCalRecHit();
  void run_kernels(const KernelConstantData<HGCeeUncalibratedRecHitConstantData>&);
  void run_kernels(const KernelConstantData<HGChefUncalibratedRecHitConstantData>&);
  void run_kernels(const KernelConstantData<HGChebUncalibratedRecHitConstantData>&);
  HGCRecHitSoA* get_output();

 private:
  void after_kernel();
  void assign_and_transfer_to_device();
  void transfer_to_host_and_synchronize();
  void reuse_device_pointers();

  const DetId::Detector dtype_;
  const std::vector<HGCUncalibratedRecHitSoA> h_oldhits_;
  KernelModifiableData<HGCUncalibratedRecHitSoA, HGCRecHitSoA> data_;
};

#endif //_KERNELMANAGER_H_
