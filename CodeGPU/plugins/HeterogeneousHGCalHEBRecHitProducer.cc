#include "HeterogeneousHGCalHEBRecHitProducer.h"

HeterogeneousHGCalHEBRecHitProducer::HeterogeneousHGCalHEBRecHitProducer(const edm::ParameterSet& ps):
  token_(consumes<HGCUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCHEBUncalibRecHitsTok")))
{
  histo1_ = fs->make<TH1F>( "energy"  , "E", 100,  0., 10. );
  histo2_ = fs->make<TH1F>( "time"  , "t", 100,  0., 10. );
  histo3_ = fs->make<TH1F>( "timeError"  , "time_error", 100,  0., 10. );
  histo4_ = fs->make<TH1I>( "son"  , "son", 100,  0., 10. );

  nhitsmax_                 = ps.getParameter<uint32_t>("nhitsmax");
  cdata_.hgcHEB_keV2DIGI_   = ps.getParameter<double>("HGCHEB_keV2DIGI");
  cdata_.hgcHEB_noise_MIP_  = ps.getParameter<edm::ParameterSet>("HGCHEB_noise_MIP").getParameter<double>("noise_MIP");
  cdata_.rangeMatch_        = ps.getParameter<uint32_t>("rangeMatch");
  cdata_.rangeMask_         = ps.getParameter<uint32_t>("rangeMask");
  cdata_.hgcHEB_isSiFE_     = ps.getParameter<bool>("HGCHEB_isSiFE");
  vdata_.weights            = ps.getParameter< std::vector<double> >("weights");
  cdata_.fhOffset_          = ps.getParameter<uint32_t>("offset"); //ddd_->layers(true);
  cdata_.s_weights_ = vdata_.weights.size();
  cdata_.hgchebUncalib2GeV_ = 1e-6 / cdata_.hgcHEB_keV2DIGI_;

  begin = std::chrono::steady_clock::now();

  tools_.reset(new hgcal::RecHitTools());
  stride_ = ( (nhitsmax_-1)/32 + 1 ) * 32; //align to warp boundary

  allocate_memory_();

  convert_constant_data_(h_kcdata_);

  produces<HGChebRecHitCollection>(collection_name_);
}

HeterogeneousHGCalHEBRecHitProducer::~HeterogeneousHGCalHEBRecHitProducer()
{
  delete kmdata_;
  delete h_kcdata_;
  delete d_kcdata_;
  delete old_soa_;
  delete d_oldhits_;
  delete d_newhits_;
  delete d_newhits_final_;
  delete h_newhits_;

  end = std::chrono::steady_clock::now();
  std::cout << "Time difference (heterogeneous) = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " [ms]" << std::endl;
}

void HeterogeneousHGCalHEBRecHitProducer::acquire(edm::Event const& event, edm::EventSetup const& setup, edm::WaitingTaskWithArenaHolder w) {
  const CUDAScopedContextAcquire ctx{event.streamID(), std::move(w), ctxState_};

  set_geometry_(setup);
  event.getByToken(token_, handle_heb_);
  const auto &hits_heb = *handle_heb_;
  unsigned int nhits = hits_heb.size();
  std::cout << "HEB hits: " << nhits << std::endl;

  convert_collection_data_to_soa_(hits_heb, old_soa_, nhits);

  kmdata_ = new KernelModifiableData<HGCUncalibratedRecHitSoA, HGCRecHitSoA>(nhitsmax_, stride_, old_soa_, d_oldhits_, d_newhits_, d_newhits_final_, h_newhits_);
  KernelManagerHGCalRecHit kernel_manager(kmdata_);

  kernel_manager.run_kernels(h_kcdata_, d_kcdata_);
  new_soa_ = kernel_manager.get_output();
  print_to_histograms(kmdata_->h_out, histo1_, histo2_, histo3_, histo4_, nhits);
  rechits_raw_ = convert_soa_data_to_collection_(kmdata_->h_out, nhits);
}

void HeterogeneousHGCalHEBRecHitProducer::produce(edm::Event& event, const edm::EventSetup& setup)
{
  CUDAScopedContextProduce ctx{ctxState_}; //only for GPU to GPU producers
  rechits_ = std::make_unique< HGCRecHitCollection >(rechits_raw_);
  event.put(std::move(rechits_), collection_name_);
}

void HeterogeneousHGCalHEBRecHitProducer::allocate_memory_()
{
  old_soa_ = new HGCUncalibratedRecHitSoA();
  d_oldhits_ = new HGCUncalibratedRecHitSoA();
  d_newhits_ = new HGCUncalibratedRecHitSoA();
  d_newhits_final_ = new HGCRecHitSoA();
  h_newhits_ = new HGCRecHitSoA();
  h_kcdata_ = new KernelConstantData<HGChebUncalibratedRecHitConstantData>(cdata_, vdata_);
  d_kcdata_ = new KernelConstantData<HGChebUncalibratedRecHitConstantData>(cdata_, vdata_);

  //_allocate pinned memory for constants on the host
  memory::allocation::host(h_kcdata_, h_mem_const_);
  //_allocate pinned memory for constants on the device
  memory::allocation::device(d_kcdata_, d_mem_const_);
  //_allocate memory for hits on the host
  memory::allocation::host(stride_, old_soa_, h_mem_in_);
  //_allocate memory for hits on the device
  memory::allocation::device(stride_, d_oldhits_, d_newhits_, d_newhits_final_, d_mem_);
  //_allocate memory for hits on the host
  memory::allocation::host(stride_, h_newhits_, h_mem_out_);
}

void HeterogeneousHGCalHEBRecHitProducer::set_geometry_(const edm::EventSetup& setup)
{
  tools_->getEventSetup(setup);
  std::string handle_str;
  handle_str = "HGCalHEScintillatorSensitive";
  edm::ESHandle<HGCalGeometry> handle;
  setup.get<IdealGeometryRecord>().get(handle_str, handle);
  ddd_ = &(handle->topology().dddConstants());
  cdata_.fhOffset_ = ddd_->layers(true);
}

void HeterogeneousHGCalHEBRecHitProducer::convert_constant_data_(KernelConstantData<HGChebUncalibratedRecHitConstantData> *kcdata)
{
  for(int i=0; i<kcdata->data.s_weights_; ++i)
    kcdata->data.weights_[i] = kcdata->vdata.weights[i];
}

void HeterogeneousHGCalHEBRecHitProducer::convert_collection_data_to_soa_(const edm::SortedCollection<HGCUncalibratedRecHit>& hits, HGCUncalibratedRecHitSoA* d, const unsigned int& nhits)
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

edm::SortedCollection<HGCRecHit> HeterogeneousHGCalHEBRecHitProducer::convert_soa_data_to_collection_(HGCRecHitSoA *d, const unsigned int& nhits)
{
  out_data_.reserve(nhits);
  for(uint i=0; i<nhits; ++i)
    {
      DetId id_converted( d->id[i] );
      out_data_[i] = HGCRecHit(id_converted, d->energy[i], d->time[i], 0, d->flagBits[i]);
    }
  return out_data_;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HeterogeneousHGCalHEBRecHitProducer);
