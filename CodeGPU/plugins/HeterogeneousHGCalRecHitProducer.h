#ifndef HeterogeneousHGCalRecHitProducer_h
#define HeterogeneousHGCalRecHitProducer_h

#include <iostream>
#include <string>
#include <memory>
#include <cuda_runtime.h>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
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

class HeterogeneousHGCalRecHitsProd: public edm::stream::EDProducer<> 
{
 public:
  explicit HeterogeneousHGCalRecHitsProd(const edm::ParameterSet& ps);
  ~HeterogeneousHGCalRecHitsProd() override;

  void produce(edm::Event& iEvent, const edm::EventSetup&) override;

 private:
  std::array< edm::EDGetTokenT<HGCUncalibratedRecHitCollection>, 3 > hits_tokens_;
  const std::array< std::string, 3 > collection_names_ = {{"HeterogeneousHGCeeUncalibratedRecHits",
 							  "HeterogeneousHGChefUncalibratedRecHits",
							  "HeterogeneousHGChebUncalibratedRecHits"}};
  edm::Handle<HGCeeUncalibratedRecHitCollection> handle_; 
};

#endif //HeterogeneousHGCalRecHitProducer_h
