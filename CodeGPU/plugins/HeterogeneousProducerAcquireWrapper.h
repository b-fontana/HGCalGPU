#ifndef _HGCALRECHITKERNELWRAPPERS_H_
#define _HGCALRECHITKERNELWRAPPERS_H_

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

#include "Utils.h"
#include "Types.h"

template <class T_IN, class T_OUT>
class HeterogeneousProducerAcquireWrapper {
 public:
  HeterogeneousProducerAcquireWrapper(const edm::SortedCollection<T_IN>&, const edm::EventSetup&);
  ~HeterogeneousProducerAcquireWrapper();
  edm::SortedCollection<T_OUT> get_output_collection();
  void run();
  void run(const CUDAScopedContextAcquire&);

 private:
  HGCUncalibratedRecHitSoA *old_soa_;
  HGCRecHitSoA *new_soa_;

  //methods
  void set_geometry_(const edm::EventSetup&);

  template <class U> void convert_collection_data_to_soa_();
  template <class U> void convert_soa_data_to_collection_();

  //geometry
  std::unique_ptr<hgcal::RecHitTools> tools_;
  const HGCalDDDConstants* ddd_;

  //data processing
  size_t size_;
  edm::SortedCollection<T_IN> hits_;
  edm::SortedCollection<T_OUT> out_data_;
  DetId::Detector det_;
};
							
#endif // _HGCALRECHITKERNELWRAPPERS_H_
