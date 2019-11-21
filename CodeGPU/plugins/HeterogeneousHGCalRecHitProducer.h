#ifndef HeterogeneousHGCalRecHitProducer_h
#define HeterogeneousHGCalRecHitProducer_h

#include <iostream>
#include <string>
#include <vector>
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

#include "HGCalRecHitKernel.cuh"

class HeterogeneousHGCalRecHitsProd: public edm::stream::EDProducer<> 
{
 public:
  explicit HeterogeneousHGCalRecHitsProd(const edm::ParameterSet& ps);
  ~HeterogeneousHGCalRecHitsProd() override;

  void produce(edm::Event& iEvent, const edm::EventSetup&) override;

 private:
  std::array< edm::EDGetTokenT<HGCUncalibratedRecHitCollection>, 3 > hitsTokens_;
  const std::array< std::string, 3 > collectionNames_ = {{"HeterogeneousHGCeeUncalibratedRecHits",
							  "HeterogeneousHGChebUncalibratedRecHits", 
							  "HeterogeneousHGChefUncalibratedRecHits"}};
  std::string message_;
  char* buffer_;
  edm::SortedCollection<HGCUncalibratedRecHit> oldhits;
  edm::SortedCollection<HGCUncalibratedRecHit> newhits;
};

#endif //HeterogeneousHGCalRecHitProducer_h
