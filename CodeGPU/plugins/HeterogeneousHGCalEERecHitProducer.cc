#include "HeterogeneousHGCalEERecHitProducer.h"
<<<<<<< HEAD
#include "HeterogeneousProducerAcquireWrapper.h"

HeterogeneousHGCalEERecHitsProd::HeterogeneousHGCalEERecHitsProd(const edm::ParameterSet& ps):
  token_(consumes<HGCUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCEEUncalibRecHitsTok")))
{
  //read all constants and prepare their move to the wrapper
  produces<HGCeeRecHitCollection>(collection_name_);
}

HeterogeneousHGCalEERecHitsProd::~HeterogeneousHGCalEERecHitsProd()
{
}

void HeterogeneousHGCalEERecHitsProd::acquire(edm::Event const& event, edm::EventSetup const& setup, edm::WaitingTaskWithArenaHolder w) {
  const CUDAScopedContextAcquire ctx{event.streamID(), std::move(w), ctxState_};
  event.getByToken(token_, handle_ee_);
  const auto &hits_ee = *handle_ee_;
  handle_size_ = hits_ee.size();

  HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit> acq_wrap(hits_ee, setup);
  acq_wrap.run();
  rechits_raw_ = acq_wrap.get_output_collection();
}

void HeterogeneousHGCalEERecHitsProd::produce(edm::Event& event, const edm::EventSetup& setup)
{
  CUDAScopedContextProduce ctx{ctxState_}; //only for GPU to GPU producers
  rechits_ = std::make_unique< HGCRecHitCollection >(rechits_raw_);
  event.put(std::move(rechits_), collection_name_);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HeterogeneousHGCalEERecHitsProd);
