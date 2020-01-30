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

#include "HeterogeneousCore/CUDACore/interface/CUDAScopedContext.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAContextState.h"
#include "HeterogeneousCore/CUDAServices/interface/CUDAService.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include "Utils.h"

class HeterogeneousHGCalEERecHitsProd: public edm::stream::EDProducer<edm::ExternalWork> 
{
 public:
  explicit HeterogeneousHGCalEERecHitsProd(const edm::ParameterSet& ps);
  ~HeterogeneousHGCalEERecHitsProd() override;

  virtual void acquire(edm::Event const&, edm::EventSetup const&, edm::WaitingTaskWithArenaHolder) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;

 private:
  edm::EDGetTokenT<HGCeeUncalibratedRecHitCollection> token_;
  const std::string collection_name_ = "HeterogeneousHGCeeUncalibratedRecHits";
  edm::Handle<HGCeeUncalibratedRecHitCollection> handle_ee_; 
  size_t handle_size_;
  HGCeeRecHitCollection rechits_raw_;
  std::unique_ptr< HGCeeRecHitCollection > rechits_;
  CUDAContextState ctxState_;

  //constants
  HGCeeUncalibratedRecHitConstantData cdata_;
};

#endif //HeterogeneousHGCalEERecHitProducer_h
