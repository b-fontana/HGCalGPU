#ifndef _HETEROGENEOUSHGCALPRODUCERACQUIREWRAPPER_H_
#define _HETEROGENEOUSHGCALPRODUCERACQUIREWRAPPER_H_

#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>
#include <type_traits>
#include <numeric>
#include <cuda_runtime.h>


#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"

#include "HeterogeneousCore/CUDACore/interface/CUDAScopedContext.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAContextState.h"
#include "HeterogeneousCore/CUDAServices/interface/CUDAService.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HeterogeneousCore/CUDAUtilities/interface/device_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/host_noncached_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/host_unique_ptr.h"

#include "KernelManager.h"
#include "Utils.h"
#include "Types.h"

template <class T_IN, class T_OUT>
class HeterogeneousHGCalProducerAcquireWrapper {
 public:
  HeterogeneousHGCalProducerAcquireWrapper(const edm::SortedCollection<T_IN>&, const edm::EventSetup&);
  ~HeterogeneousHGCalProducerAcquireWrapper();
  edm::SortedCollection<T_OUT> get_output_collection();
  template <class U> void run(KernelConstantData<U>&);
  //template <class U> void run(const KernelConstantData<U>&, const CUDAScopedContextAcquire&);

 private:
  HGCUncalibratedRecHitSoA *old_soa_ = nullptr, *d_oldhits_ = nullptr, *d_newhits_ = nullptr;
  HGCRecHitSoA *new_soa_ = nullptr, *d_newhits_final_ = nullptr, *h_newhits_ = nullptr;

  //methods
  void set_geometry_(const edm::EventSetup&);
  std::tuple<size_t, size_t, size_t> get_memory_sizes_(const std::vector<size_t>&, const size_t&, const size_t&);
  void allocate_host_(KernelConstantData<HGCeeUncalibratedRecHitConstantData>&, cudautils::host::noncached::unique_ptr<double[]>&);
  void allocate_host_(KernelConstantData<HGChefUncalibratedRecHitConstantData>&, cudautils::host::noncached::unique_ptr<double[]>&);
  void allocate_host_(KernelConstantData<HGChebUncalibratedRecHitConstantData>&, cudautils::host::noncached::unique_ptr<double[]>&);
  void allocate_host_(HGCUncalibratedRecHitSoA*&, cudautils::host::noncached::unique_ptr<float[]>&);
  void allocate_host_(HGCRecHitSoA*&, cudautils::host::unique_ptr<float[]>&);
  void allocate_device_(KernelConstantData<HGCeeUncalibratedRecHitConstantData>&, cudautils::device::unique_ptr<double[]>&);
  void allocate_device_(KernelConstantData<HGChefUncalibratedRecHitConstantData>&, cudautils::device::unique_ptr<double[]>&);
  void allocate_device_(KernelConstantData<HGChebUncalibratedRecHitConstantData>&, cudautils::device::unique_ptr<double[]>&);
  void allocate_device_(HGCUncalibratedRecHitSoA*&, HGCUncalibratedRecHitSoA*&, HGCRecHitSoA*&, cudautils::device::unique_ptr<float[]>&);
  template <class U> void convert_collection_data_to_soa_();
  template <class U> void convert_soa_data_to_collection_();
  void convert_constant_data_(KernelConstantData<HGCeeUncalibratedRecHitConstantData>&);
  void convert_constant_data_(KernelConstantData<HGChefUncalibratedRecHitConstantData>&);
  void convert_constant_data_(KernelConstantData<HGChebUncalibratedRecHitConstantData>&);

  //geometry
  std::unique_ptr<hgcal::RecHitTools> tools_;
  const HGCalDDDConstants* ddd_ = nullptr;

  //data processing
  size_t nhits_;
  unsigned int stride_;
  edm::SortedCollection<T_IN> hits_;
  DetId::Detector dtype_;
  edm::SortedCollection<T_OUT> out_data_;
  DetId::Detector det_;
};
							
#endif // _HETEROGENEOUSHGCALPRODUCERACQUIREWRAPPER_H_
