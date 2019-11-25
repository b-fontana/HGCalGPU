#include "HeterogeneousHGCalRecHitProducer.h"

HeterogeneousHGCalRecHitsProd::HeterogeneousHGCalRecHitsProd(const edm::ParameterSet& ps):
  hits_tokens_({{consumes<HGCeeUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCEEUncalibRecHitsTok")),
	  consumes<HGChefUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCHEFUncalibRecHitsTok")),
	  consumes<HGChebUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCHEBUncalibRecHitsTok"))}})
{
  for(unsigned int idet=Detector::EE; idet<Detector::NDetectors; ++idet)
      produces<HGCUncalibratedRecHitCollection>(collection_names_[idet]);
}

HeterogeneousHGCalRecHitsProd::~HeterogeneousHGCalRecHitsProd()
{
}

void HeterogeneousHGCalRecHitsProd::produce(edm::Event& iEvent, const edm::EventSetup&)
{
  for(unsigned int idet=Detector::EE; idet<Detector::NDetectors; ++idet)
    {
      iEvent.getByToken(hits_tokens_[idet], handle_);
      const auto &hits = *handle_;
      auto rechits = hgc_rechit_kernel_wrapper<HGCUncalibratedRecHitCollection>(hits, static_cast<Detector>(idet));
      for(const auto &hit: hits)
	rechits->push_back(hit);  
      iEvent.put(std::move(rechits), collection_names_[idet]);
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HeterogeneousHGCalRecHitsProd);
