#include "HeterogeneousHGCalEERecHitProducer.h"
#include "HeterogeneousHGCalProducerAcquireWrapper.h"

HeterogeneousHGCalEERecHitsProd::HeterogeneousHGCalEERecHitsProd(const edm::ParameterSet& ps):
  token_(consumes<HGCUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCEEUncalibRecHitsTok")))
{
  cdata_.hgcEE_keV2DIGI_   = ps.getParameter<double>("HGCEE_keV2DIGI");
  cdata_.rangeMatch_       = ps.getParameter<uint32_t>("rangeMatch");
  cdata_.rangeMask_        = ps.getParameter<uint32_t>("rangeMask");
  cdata_.hgcEE_isSiFE_     = ps.getParameter<bool>("HGCEE_isSiFE");
  vdata_.fCPerMIP          = ps.getParameter< std::vector<double> >("HGCEE_fCPerMIP");
  vdata_.cce               = ps.getParameter<edm::ParameterSet>("HGCEE_cce").getParameter<std::vector<double> >("values");
  vdata_.noise_fC          = ps.getParameter<edm::ParameterSet>("HGCEE_noise_fC").getParameter<std::vector<double> >("values");
  vdata_.rcorr             = ps.getParameter< std::vector<double> >("rcorr");
  vdata_.weights           = ps.getParameter< std::vector<double> >("weights");

  cdata_.s_hgcEE_fCPerMIP_ = vdata_.fCPerMIP.size();
  cdata_.s_hgcEE_cce_      = vdata_.cce.size();
  cdata_.s_hgcEE_noise_fC_ = vdata_.noise_fC.size();
  cdata_.s_rcorr_            = vdata_.rcorr.size();
  cdata_.s_weights_ = vdata_.weights.size();

  cdata_.hgceeUncalib2GeV_ = 1e-6 / cdata_.hgcEE_keV2DIGI_;
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
