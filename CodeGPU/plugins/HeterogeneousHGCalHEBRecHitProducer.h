#ifndef HeterogeneousHGCalHEBRecHitProducer_h
#define HeterogeneousHGCalHEBRecHitProducer_h

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

#include "HGCalRecHitKernelWrappers.h"
#include "Utils.h"

class HeterogeneousHGCalHEBRecHitsProd: public edm::stream::EDProducer<> 
{
 public:
  explicit HeterogeneousHGCalHEBRecHitsProd(const edm::ParameterSet& ps);
  ~HeterogeneousHGCalHEBRecHitsProd() override;

  void produce(edm::Event& iEvent, const edm::EventSetup&) override;

 private:
  edm::EDGetTokenT<HGChebUncalibratedRecHitCollection> token_;
  const std::string collection_name_ = "HeterogeneousHGChebUncalibratedRecHits";
  edm::Handle<HGChebUncalibratedRecHitCollection> handle_heb_;  
};

#endif //HeterogeneousHGCalHEBRecHitProducer_h
