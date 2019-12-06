#ifndef HeterogeneousHGCalHEFRecHitProducer_h
#define HeterogeneousHGCalHEFRecHitProducer_h

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

class HeterogeneousHGCalHEFRecHitsProd: public edm::stream::EDProducer<> 
{
 public:
  explicit HeterogeneousHGCalHEFRecHitsProd(const edm::ParameterSet& ps);
  ~HeterogeneousHGCalHEFRecHitsProd() override;

  void produce(edm::Event& iEvent, const edm::EventSetup&) override;

 private:
  edm::EDGetTokenT<HGChefUncalibratedRecHitCollection> token_;
  const std::string collection_name_ = "HeterogeneousHGChefUncalibratedRecHits";
  edm::Handle<HGChefUncalibratedRecHitCollection> handle_hef_;
};

#endif //HeterogeneousHGCalHEFRecHitProducer_h
