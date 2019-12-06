#ifndef _KERNELMANAGER_H
#define _KERNELMANAGER_H

#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
//#include "DataFormats/​HGCRecHit/​interface/​HGCRecHit.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCompat.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "Utils.h"

#include <algorithm> //std::swap  
#include <variant>
#include <cuda.h>
#include <cuda_runtime.h>

#ifdef __CUDA_ARCH__
extern __constant__ uint32_t calo_rechit_masks[];
#endif

class KernelManagerBase  {
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
  explicit KernelManagerHGCalRecHit(const edm::SortedCollection<HGCUncalibratedRecHit>&);
  ~KernelManagerHGCalRecHit();
  void run_kernels();
  edm::SortedCollection<HGCRecHit> get_new_collection();

private:
  void ee_step1_wrapper();
  void hef_step1_wrapper();
  void heb_step1_wrapper();
  void to_rechit_wrapper();
  void assign_and_transfer_to_device() override;
  void transfer_to_host_and_synchronize() override;
  void reuse_device_pointers() override;

  size_t shits_, sbytes_;
  DetId::Detector dtype_;
  edm::SortedCollection<HGCUncalibratedRecHit> oldhits_collection_;
  HGCUncalibratedRecHit *h_oldhits_; //host pointer
  HGCRecHit *h_newhits_; //host pointer
  HGCUncalibratedRecHit *d_oldhits_, *d_newhits_; //device pointers
  HGCRecHit *d_newhits_final_; //device pointer
};

#endif //_KERNELMANAGER_H_
