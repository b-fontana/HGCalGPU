#include "HeterogeneousHGCalRecHitProducer.h"

HeterogeneousHGCalRecHitsProd::HeterogeneousHGCalRecHitsProd(const edm::ParameterSet& ps):
  hitsTokens_({{consumes<HGCeeUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCEEUncalibRecHitsTok")),
	  consumes<HGChefUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCHEFUncalibRecHitsTok")),
	  consumes<HGChebUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCHEBUncalibRecHitsTok"))}}),
  message_(ps.getUntrackedParameter<std::string>("message"))
{
  produces<HGCeeUncalibratedRecHitCollection>(collectionNames_[0]);
  produces<HGCeeUncalibratedRecHitCollection>(collectionNames_[1]);
  produces<HGCeeUncalibratedRecHitCollection>(collectionNames_[2]);
  cudaCheck(cudaMalloc(&buffer_, message_.size()));
  cudaCheck(cudaMemcpy(buffer_, message_.data(), message_.size(), cudaMemcpyDefault));
}

HeterogeneousHGCalRecHitsProd::~HeterogeneousHGCalRecHitsProd()
{
  cudaCheck(cudaFree(buffer_));
}

void HeterogeneousHGCalRecHitsProd::produce(edm::Event& iEvent, const edm::EventSetup&)
{
  scrambler_wrapper(buffer_, message_.size());
  
  auto hefRecHits = std::make_unique<HGChefUncalibratedRecHitCollection>();
  auto hebRecHits = std::make_unique<HGChebUncalibratedRecHitCollection>();

  edm::Handle<HGCeeUncalibratedRecHitCollection> eeRecHitsHandle;
  edm::Handle<HGChefUncalibratedRecHitCollection> hefRecHitsHandle;
  edm::Handle<HGChebUncalibratedRecHitCollection> hebRecHitsHandle;

  iEvent.getByToken(hitsTokens_[0], eeRecHitsHandle);
  const auto &eeHits = *eeRecHitsHandle;
  std::cout << "check" << std::endl;
  auto eeRecHits = HGCeeRecHitKernel_wrapper(eeHits, collectionNames_[0]);
  std::cout << "check" << std::endl;

  iEvent.getByToken(hitsTokens_[1], hefRecHitsHandle);
  const auto &hefHits = *hefRecHitsHandle;  
  for(const auto &hit: hefHits)
      hefRecHits->push_back(hit);

  iEvent.getByToken(hitsTokens_[2], hebRecHitsHandle);
  const auto &hebHits = *hebRecHitsHandle;  
  for(const auto &hit: hebHits)
      hebRecHits->push_back(hit);

  iEvent.put(std::move(eeRecHits), collectionNames_[0]); //product instance label
  iEvent.put(std::move(hefRecHits), collectionNames_[1]); //product instance label
  iEvent.put(std::move(hebRecHits), collectionNames_[2]); //product instance label
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HeterogeneousHGCalRecHitsProd);
