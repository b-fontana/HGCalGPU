#include "HeterogeneousHGCalEERecHitProducer.h"
#include "HGCalRecHitKernelWrappers.h"

HeterogeneousHGCalEERecHitsProd::HeterogeneousHGCalEERecHitsProd(const edm::ParameterSet& ps):
  token_(consumes<HGCUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCEEUncalibRecHitsTok")))
{
  //read all constants and prepare their move to the wrapper
  produces<HGCeeUncalibratedRecHitCollection>(collection_name_);
}

HeterogeneousHGCalEERecHitsProd::~HeterogeneousHGCalEERecHitsProd()
{
}

void HeterogeneousHGCalEERecHitsProd::acquire(edm::Event const& event, edm::EventSetup const& setup, edm::WaitingTaskWithArenaHolder w) {
  const CUDAScopedContextAcquire ctx{event.streamID(), std::move(w), ctxState_};

  event.getByToken(token_, handle_ee_);
  const auto &hits_ee = *handle_ee_;
  handle_size_ = hits_ee.size();

  //call 'set()' here through a templated 'geometry_wrapper' stored next to kernel_manager_wrapper()
  //the same wrapper should work for all producers, ideally even for the Digi phase

  rechits_raw_ = kernel_manager_wrapper<HGCUncalibratedRecHit, HGCRecHit>(hits_ee, ctx);
}

void HeterogeneousHGCalEERecHitsProd::produce(edm::Event& event, const edm::EventSetup& setup)
{
  CUDAScopedContextProduce ctx{ctxState_};

  rechits_ = std::make_unique< HGCRecHitCollection >(rechits_raw_);
  event.put(std::move(rechits_), collection_name_);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HeterogeneousHGCalEERecHitsProd);
