#ifndef HeterogeneousHGCalEERecHitProducer_h
#define HeterogeneousHGCalEERecHitProducer_h

#include <iostream>
#include <string>
#include <memory>
#include <cuda_runtime.h>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"

#include "HeterogeneousCore/CUDACore/interface/CUDAScopedContext.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAContextState.h"
#include "HeterogeneousCore/CUDAServices/interface/CUDAService.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include "HeterogeneousHGCalProducerMemoryWrapper.h"
#include "KernelManager.h"
#include "Utils.h"

class HeterogeneousHGCalEERecHitsProd: public edm::stream::EDProducer<edm::ExternalWork> 
{
 public:
  explicit HeterogeneousHGCalEERecHitsProd(const edm::ParameterSet& ps);
  ~HeterogeneousHGCalEERecHitsProd() override;

  virtual void acquire(edm::Event const&, edm::EventSetup const&, edm::WaitingTaskWithArenaHolder) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;

 private:
  unsigned int nhitsmax_ = 0;
  unsigned int stride_ = 0;
  edm::EDGetTokenT<HGCeeUncalibratedRecHitCollection> token_;
  const std::string collection_name_ = "HeterogeneousHGCeeUncalibratedRecHits";
  edm::Handle<HGCeeUncalibratedRecHitCollection> handle_ee_; 
  size_t handle_size_;
  HGCeeRecHitCollection rechits_raw_;
  std::unique_ptr< HGCeeRecHitCollection > rechits_;
  CUDAContextState ctxState_;

  //constants
  HGCeeUncalibratedRecHitConstantData cdata_;
  HGCConstantVectorData vdata_;

  //memory
  cudautils::host::noncached::unique_ptr<double[]> h_mem_const_;
  cudautils::device::unique_ptr<double[]> d_mem_const_;
  cudautils::host::noncached::unique_ptr<float[]> h_mem_in_;
  cudautils::device::unique_ptr<float[]> d_mem_;
  cudautils::host::unique_ptr<float[]> h_mem_out_;

  //geometry
  void set_geometry_(const edm::EventSetup&);
  std::unique_ptr<hgcal::RecHitTools> tools_;
  const HGCalDDDConstants* ddd_ = nullptr;

  //data processing
  void convert_collection_data_to_soa_(const edm::SortedCollection<HGCUncalibratedRecHit>&, HGCUncalibratedRecHitSoA*, const unsigned int&);
  edm::SortedCollection<HGCRecHit> convert_soa_data_to_collection_(HGCRecHitSoA*, const unsigned int&);
  void convert_constant_data_(KernelConstantData<HGCeeUncalibratedRecHitConstantData>*);

  HGCUncalibratedRecHitSoA *old_soa_ = nullptr, *d_oldhits_ = nullptr, *d_newhits_ = nullptr;
  HGCRecHitSoA *new_soa_ = nullptr, *d_newhits_final_ = nullptr, *h_newhits_ = nullptr;
  KernelModifiableData<HGCUncalibratedRecHitSoA, HGCRecHitSoA> *kmdata_;
  KernelConstantData<HGCeeUncalibratedRecHitConstantData> *h_kcdata_;
  KernelConstantData<HGCeeUncalibratedRecHitConstantData> *d_kcdata_;
  edm::SortedCollection<HGCRecHit> out_data_;

};

#endif //HeterogeneousHGCalEERecHitProducer_h
