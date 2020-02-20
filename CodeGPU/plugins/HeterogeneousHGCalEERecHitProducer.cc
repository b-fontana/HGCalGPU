#include "HeterogeneousHGCalEERecHitProducer.h"

HeterogeneousHGCalEERecHitsProd::HeterogeneousHGCalEERecHitsProd(const edm::ParameterSet& ps):
  token_(consumes<HGCUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCEEUncalibRecHitsTok")))
{
  nhitsmax_                = ps.getParameter<uint32_t>("nhitsmax");
  cdata_.hgcEE_keV2DIGI_   = ps.getParameter<double>("HGCEE_keV2DIGI");
  cdata_.xmin_             = ps.getParameter<double>("minValSiPar"); //float
  cdata_.xmax_             = ps.getParameter<double>("maxValSiPar"); //float
  cdata_.aterm_            = ps.getParameter<double>("constSiPar"); //float
  cdata_.cterm_            = ps.getParameter<double>("noiseSiPar"); //float
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
  vdata_.waferTypeL = {0, 1, 2};//ddd_->retWaferTypeL(); if depends on geometry the allocation is tricky!
  cdata_.s_waferTypeL_ = vdata_.waferTypeL.size();
  //cdata.data.fhOffset = ddd_->layers(true);

  tools_.reset(new hgcal::RecHitTools());
  stride_ = ( (nhitsmax_-1)/32 + 1 ) * 32; //align to warp boundary
  old_soa_ = new HGCUncalibratedRecHitSoA();
  d_oldhits_ = new HGCUncalibratedRecHitSoA();
  d_newhits_ = new HGCUncalibratedRecHitSoA();
  d_newhits_final_ = new HGCRecHitSoA();
  h_newhits_ = new HGCRecHitSoA();
  h_kcdata_ = new KernelConstantData<HGCeeUncalibratedRecHitConstantData>(cdata_, vdata_);
  std::cout << "Address outside before: " << d_oldhits_->amplitude << std::endl;
  d_kcdata_ = new KernelConstantData<HGCeeUncalibratedRecHitConstantData>(cdata_, vdata_);
  std::cout << "Address outside before: " << old_soa_->amplitude << std::endl;
  std::cout << "Address outside before: " << d_oldhits_->amplitude << std::endl;

  //_allocate pinned memory for constants on the host
  memory::allocation::host(h_kcdata_, h_mem_const_);
  //_allocate pinned memory for constants on the device
  memory::allocation::device(d_kcdata_, d_mem_const_);
  //_allocate memory for hits on the host
  memory::allocation::host(nhitsmax_, old_soa_, h_mem_in_);
  //_allocate memory for hits on the device
  memory::allocation::device(nhitsmax_, d_oldhits_, d_newhits_, d_newhits_final_, d_mem_);
  //_allocate memory for hits on the host
  memory::allocation::host(nhitsmax_, h_newhits_, h_mem_out_);

  std::cout << "Address outside after: " << old_soa_->amplitude << std::endl;
  std::cout << "Address outside after: " << d_oldhits_->amplitude << std::endl;

  convert_constant_data_(h_kcdata_);
  produces<HGCeeRecHitCollection>(collection_name_);
}

HeterogeneousHGCalEERecHitsProd::~HeterogeneousHGCalEERecHitsProd()
{
  delete kmdata_;
  delete h_kcdata_;
  delete d_kcdata_;
  delete old_soa_;
  delete d_oldhits_;
  delete d_newhits_;
  delete d_newhits_final_;
  delete h_newhits_;
}

void HeterogeneousHGCalEERecHitsProd::acquire(edm::Event const& event, edm::EventSetup const& setup, edm::WaitingTaskWithArenaHolder w) {
  const CUDAScopedContextAcquire ctx{event.streamID(), std::move(w), ctxState_};

  set_geometry_(setup);
  event.getByToken(token_, handle_ee_);
  const auto &hits_ee = *handle_ee_;
  unsigned int nhits = hits_ee.size();
  convert_collection_data_to_soa_(hits_ee, old_soa_, nhits);

  kmdata_ = new KernelModifiableData<HGCUncalibratedRecHitSoA, HGCRecHitSoA>(nhitsmax_, stride_, old_soa_, d_oldhits_, d_newhits_, d_newhits_final_, h_newhits_);
  KernelManagerHGCalRecHit kernel_manager(kmdata_);

  kernel_manager.run_kernels(h_kcdata_, d_kcdata_);
  std::cout << "check3"  << std::endl;
  new_soa_ = kernel_manager.get_output();
  std::cout << "check4"  << std::endl;
  rechits_raw_ = convert_soa_data_to_collection_(kmdata_->h_out, nhits);
  std::cout << "check5"  << std::endl;
}

void HeterogeneousHGCalEERecHitsProd::produce(edm::Event& event, const edm::EventSetup& setup)
{
  CUDAScopedContextProduce ctx{ctxState_}; //only for GPU to GPU producers
  std::cout << "check6"  << std::endl;
  rechits_ = std::make_unique< HGCRecHitCollection >(rechits_raw_);
  std::cout << "check7"  << std::endl;
  event.put(std::move(rechits_), collection_name_);
}

void HeterogeneousHGCalEERecHitsProd::set_geometry_(const edm::EventSetup& setup)
{
  tools_->getEventSetup(setup);
  //rechitMaker_->set(es);
  std::string handle_str;
  handle_str = "HGCalEESensitive";
  /*
  else if (det_ == DetId::HGCalHSi) 
    handle_str = "HGCalHESiliconSensitive";
  else if (det_ == DetId::HGCalHSc)
    handle_str = "HGCalHESiliconSensitive";
  */
  edm::ESHandle<HGCalGeometry> handle;
  setup.get<IdealGeometryRecord>().get(handle_str, handle);
  ddd_ = &(handle->topology().dddConstants());
}

void HeterogeneousHGCalEERecHitsProd::convert_constant_data_(KernelConstantData<HGCeeUncalibratedRecHitConstantData> *kcdata)
{
  for(int i=0; i<kcdata->data.s_hgcEE_fCPerMIP_; ++i)
    kcdata->data.hgcEE_fCPerMIP_[i] = kcdata->vdata.fCPerMIP[i];
  for(int i=0; i<kcdata->data.s_hgcEE_cce_; ++i)
    kcdata->data.hgcEE_cce_[i] = kcdata->vdata.cce[i];
  for(int i=0; i<kcdata->data.s_hgcEE_noise_fC_; ++i)
    kcdata->data.hgcEE_noise_fC_[i] = kcdata->vdata.noise_fC[i];
  for(int i=0; i<kcdata->data.s_rcorr_; ++i)
    kcdata->data.rcorr_[i] = kcdata->vdata.rcorr[i];
  for(int i=0; i<kcdata->data.s_weights_; ++i)
    kcdata->data.weights_[i] = kcdata->vdata.weights[i];
  for(int i=0; i<kcdata->data.s_waferTypeL_; ++i)
    kcdata->data.waferTypeL_[i] = kcdata->vdata.waferTypeL[i];
}

void HeterogeneousHGCalEERecHitsProd::convert_collection_data_to_soa_(const edm::SortedCollection<HGCUncalibratedRecHit>& hits, HGCUncalibratedRecHitSoA* d, const unsigned int& nhits)
{
  for(unsigned int i=0; i<nhits; ++i)
    {
      d->amplitude[i] = hits[i].amplitude();
      d->pedestal[i] = hits[i].pedestal();
      d->jitter[i] = hits[i].jitter();
      d->chi2[i] = hits[i].chi2();
      d->OOTamplitude[i] = hits[i].outOfTimeEnergy();
      d->OOTchi2[i] = hits[i].outOfTimeChi2();
      d->flags[i] = hits[i].flags();
      d->aux[i] = 0;
      d->id[i] = hits[i].id().rawId();
    }
}

edm::SortedCollection<HGCRecHit> HeterogeneousHGCalEERecHitsProd::convert_soa_data_to_collection_(HGCRecHitSoA *d, const unsigned int& nhits)
{
  out_data_.reserve(nhits);
  for(uint i=0; i<nhits; ++i)
    {
      DetId id_converted( d->id[i] );
      std::cout << "NEWSOA: " << d->id[i] << ", " << d->id[i] << ", " << d->energy[i] << ", " << d->time[i] << ", " << d->flagBits[i] << std::endl;
      out_data_[i] = HGCRecHit(id_converted, d->energy[i], d->time[i], 0/*d->flags[i]*/, d->flagBits[i]);
    }
  return out_data_;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HeterogeneousHGCalEERecHitsProd);
