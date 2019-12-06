#include "HeterogeneousHGCalHEBRecHitProducer.h"

HeterogeneousHGCalHEBRecHitsProd::HeterogeneousHGCalHEBRecHitsProd(const edm::ParameterSet& ps):
  token_(consumes<HGChebUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCHEBUncalibRecHitsTok")))
{
  produces<HGChebUncalibratedRecHitCollection>(collection_name_);
}

HeterogeneousHGCalHEBRecHitsProd::~HeterogeneousHGCalHEBRecHitsProd()
{
}

void HeterogeneousHGCalHEBRecHitsProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  iEvent.getByToken(token_, handle_heb_);
  const auto &hits_heb = *handle_heb_;

  auto rechits = kernel_manager_wrapper<HGCUncalibratedRecHit, HGCRecHit>(hits_heb, iSetup);

  iEvent.put(std::move(rechits), collection_name_);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HeterogeneousHGCalHEBRecHitsProd);
