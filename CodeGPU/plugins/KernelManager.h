#ifndef _KERNELMANAGER_H
#define _KERNELMANAGER_H

//#include "FWCore/Utilities/interface/Exception.h"
//#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
//#include "DataFormats/DetId/interface/DetId.h"
//#include "DataFormats/HGCRecHit_GPU/interface/HGCRecHitCollections.h"
//#include "DataFormats/​HGCRecHit_GPU/​interface/​HGCRecHit.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCompat.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
//#include "HeterogeneousCore/CUDACore/interface/CUDAScopedContext.h"
//#include "HeterogeneousCore/CUDAUtilities/interface/host_noncached_unique_ptr.h"
//#include "HeterogeneousCore/CUDAUtilities/interface/device_unique_ptr.h"
//#include "HeterogeneousCore/CUDAUtilities/interface/host_unique_ptr.h"
//#include "Utils.h"

#include "HGCalRecHitKernelImpl.cuh"

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
  class KernelManagerData {
 public:
  KernelManagerData(TYPE_IN *h_in, TYPE_IN *d_1, TYPE_IN *d_2,
		    TYPE_OUT *d_out, TYPE_OUT *h_out):
  h_in_(h_in), d_1_(d_1), d_2_(d_2), d_out_(d_out), h_out_(h_out) {}

  TYPE_IN *h_in_;
  TYPE_IN *d_1_, *d_2_;
  TYPE_OUT *d_out_;
  TYPE_OUT *h_out_;
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
  explicit KernelManagerHGCalRecHit(KernelManagerData<HGCUncalibratedRecHit_GPU, HGCRecHit_GPU>);
  ~KernelManagerHGCalRecHit();
  void run_kernels();
  HGCRecHit_GPU* get_output();

 private:
  friend class KernelManagerData<HGCUncalibratedRecHit_GPU, HGCRecHit_GPU>;

  void ee_step1_wrapper();
  void hef_step1_wrapper();
  void heb_step1_wrapper();
  void to_rechit_wrapper();
  void assign_and_transfer_to_device() override;
  void transfer_to_host_and_synchronize() override;
  void reuse_device_pointers() override;

  size_t shits_, sbytes_;
  //DetId::Detector dtype_;
  const std::vector<HGCUncalibratedRecHit_GPU> h_oldhits_;
  KernelManagerData<HGCUncalibratedRecHit_GPU, HGCRecHit_GPU> data_;
};

#endif //_KERNELMANAGER_H_
