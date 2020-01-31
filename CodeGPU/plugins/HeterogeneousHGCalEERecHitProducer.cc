#include "HeterogeneousHGCalEERecHitProducer.h"
#include "HeterogeneousHGCalProducerAcquireWrapper.h"

HeterogeneousHGCalEERecHitsProd::HeterogeneousHGCalEERecHitsProd(const edm::ParameterSet& ps):
  token_(consumes<HGCUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCEEUncalibRecHitsTok")))
{
  cdata_.hgcEE_keV2DIGI_   = ps.getParameter<double>("HGCEE_keV2DIGI");
  cdata_.hgceeUncalib2GeV_ = ps.getParameter<double>("   ");
  cdata_.rangeMatch_       = ps.getParameter<uint32_t>("   ");
  cdata_.rangeMask_        = ps.getParameter<uint32_t>("   ");
  cdata_.hgcEE_isSiFE_     = ps.getParameter<bool>("   ");
  vdata_.fCPerMIP      = ps.getParameter<>("   ");
  vdata_.cee           = ps.getParameter<>("   ");
  vdata_.noise_fC      = ps.getParameter<>("   ");
  vdata_.rcorr   = ps.getParameter<std::vector<float>>("   ");
  vdata_.weights = ps.getParameter<std::vector<float>>("   ");
  vdata_.sizes[0] = vdata_.fCPerMIP.size();
  vdata_.sizes[1] = vdata_.cee.size();
  vdata_.sizes[2] = vdata_.noise_fC.size();
  vdata_.sizes[3] = vdata_.rcorr.size();
  vdata_.sizes[4] = vdata_.weights.size();
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

  KernelConstantData<HGCeeUncalibratedRecHitConstantData> kcdata(cdata_, vdata_);
  HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit> acq_wrap(hits_ee, setup);
  acq_wrap.run(kcdata);
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
