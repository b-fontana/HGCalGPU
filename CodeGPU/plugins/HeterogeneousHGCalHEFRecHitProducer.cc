#include "HeterogeneousHGCalHEFRecHitProducer.h"

HeterogeneousHGCalHEFRecHitsProd::HeterogeneousHGCalHEFRecHitsProd(const edm::ParameterSet& ps):
  token_(consumes<HGChefUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCHEFUncalibRecHitsTok")))
{
  produces<HGChefUncalibratedRecHitCollection>(collection_name_);
}

HeterogeneousHGCalHEFRecHitsProd::~HeterogeneousHGCalHEFRecHitsProd()
{
}

void HeterogeneousHGCalHEFRecHitsProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  iEvent.getByToken(token_, handle_hef_);
  const auto &hits_hef = *handle_hef_;

  auto rechits = kernel_manager_wrapper<HGCUncalibratedRecHit, HGCRecHit>(hits_hef, iSetup);

  iEvent.put(std::move(rechits), collection_name_);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HeterogeneousHGCalHEFRecHitsProd);
