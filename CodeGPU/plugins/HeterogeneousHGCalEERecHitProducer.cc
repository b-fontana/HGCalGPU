#include "HeterogeneousHGCalEERecHitProducer.h"
#include "HGCalRecHitKernelWrappers.h"

HeterogeneousHGCalEERecHitsProd::HeterogeneousHGCalEERecHitsProd(const edm::ParameterSet& ps):
  token_(consumes<HGCeeUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCEEUncalibRecHitsTok")))
{
  produces<HGCeeUncalibratedRecHitCollection>(collection_name_);
}

HeterogeneousHGCalEERecHitsProd::~HeterogeneousHGCalEERecHitsProd()
{
}

void HeterogeneousHGCalEERecHitsProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  iEvent.getByToken(token_, handle_ee_);
  const auto &hits_ee = *handle_ee_;

  //call 'set()' here through a templated 'geometry_wrapper' stored next to kernel_manager_wrapper()
  //the same wrapper should work for all producers, ideally even for the Digi phase

  auto rechits = kernel_manager_wrapper<HGCUncalibratedRecHit, HGCRecHit>(hits_ee, iSetup);
      
  iEvent.put(std::move(rechits), collection_name_);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HeterogeneousHGCalEERecHitsProd);
