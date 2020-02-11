#include "HeterogeneousHGCalProducerAcquireWrapper.h"
#include "Types.h"

template <class T_IN, class T_OUT>
HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::HeterogeneousHGCalProducerAcquireWrapper(const edm::SortedCollection<T_IN>& hits, const edm::EventSetup& setup)
{
  nhits_ = hits.size();
  if (nhits_ == 0)
    throw cms::Exception("EmptyCollection") << "The passed collection is empty."; 
 
  for(unsigned int i=0; i<nhits_-1; ++i)
    assert(hits[i].id().det() == hits[i+1].id().det());
  
  det_ = hits[0].id().det(); 
  if( !(det_ == DetId::HGCalEE or det_ == DetId::HGCalHSi or det_ == DetId::HGCalHSc) )
    throw cms::Exception("WrongDetectorType") << "The specified detector is wrong.";

  stride_ = ( (nhits_-1)/32 + 1 ) * 32;
  hits_ = hits;
  out_data_.reserve(nhits_);
  tools_.reset(new hgcal::RecHitTools());
  this->set_geometry_(setup);
}

template <class T_IN, class T_OUT>
HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::~HeterogeneousHGCalProducerAcquireWrapper()
{
  delete d_oldhits_;
  delete d_newhits_;
  delete d_newhits_final_;
  delete h_newhits_;
}

template <class T_IN, class T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::set_geometry_(const edm::EventSetup& setup)
{
  tools_->getEventSetup(setup);
  //rechitMaker_->set(es);
  std::string handle_str;
  if (det_ == DetId::HGCalEE) 
    handle_str = "HGCalEESensitive";
  else if (det_ == DetId::HGCalHSi) 
    handle_str = "HGCalHESiliconSensitive";
  else if (det_ == DetId::HGCalHSc)
    handle_str = "HGCalHESiliconSensitive";
  edm::ESHandle<HGCalGeometry> handle;
  setup.get<IdealGeometryRecord>().get(handle_str, handle);
  ddd_ = &(handle->topology().dddConstants());
}

template <class T_IN, class T_OUT>
template <class U>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::set_geometry_constants_(KernelConstantData<U>& kcdata)
{
  kcdata.vdata.waferTypeL = ddd_->retWaferTypeL();
  kcdata.data.s_waferTypeL_ = kcdata.vdata.waferTypeL.size();
}

//returns total number of bytes, number of 'double' elements and number of 'float' elements
template <typename T_IN, typename T_OUT>
std::tuple<size_t, size_t, size_t, size_t> HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::get_memory_sizes_(const std::vector<size_t>& fixed_sizes, const size_t& ndoubles, const size_t& nfloats, const size_t& nints)
{
  const size_t size1 = sizeof(double);
  const size_t size2 = sizeof(float);
  const size_t size3 = sizeof(int);
  size_t nelements1_tot = std::accumulate( fixed_sizes.begin(), fixed_sizes.begin() + ndoubles, 0);
  size_t nelements2_tot = std::accumulate( fixed_sizes.begin() + ndoubles, fixed_sizes.begin() + ndoubles + nfloats, 0);
  size_t nelements3_tot = std::accumulate( fixed_sizes.begin() + ndoubles + nfloats, fixed_sizes.end(), 0);
  assert( fixed_sizes.begin() + ndoubles + nfloats + nints == fixed_sizes.end() );
  size_t size_tot = nelements1_tot*size1+nelements2_tot*size2+nelements3_tot*size3;
  return std::make_tuple(size_tot, nelements1_tot, nelements2_tot, nelements3_tot);
}

template <typename T_IN, typename T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_device_(KernelConstantData<HGCeeUncalibratedRecHitConstantData>& kcdata, cudautils::device::unique_ptr<double[]>& mem) {
  const std::vector<size_t> nelements = {kcdata.data.s_hgcEE_fCPerMIP_, kcdata.data.s_hgcEE_cce_, kcdata.data.s_hgcEE_noise_fC_, kcdata.data.s_rcorr_, kcdata.data.s_weights_, kcdata.data.s_waferTypeL_};
  auto memsizes = get_memory_sizes_(nelements, 5, 0, 1);

  mem = cudautils::make_device_unique<double[]>(std::get<0>(memsizes), 0);

  kcdata.data.hgcEE_fCPerMIP_ = mem.get();
  kcdata.data.hgcEE_cce_      = kcdata.data.hgcEE_fCPerMIP_ + nelements[0];
  kcdata.data.hgcEE_noise_fC_ = kcdata.data.hgcEE_cce_ + nelements[1];
  kcdata.data.rcorr_          = kcdata.data.hgcEE_noise_fC_ + nelements[2];
  kcdata.data.weights_        = kcdata.data.rcorr_ + nelements[3];
  kcdata.data.waferTypeL_     = reinterpret_cast<int*>(kcdata.data.weights_ + nelements[4]);
  kcdata.data.nbytes = std::get<0>(memsizes);
  kcdata.data.ndelem = std::get<1>(memsizes) + 2;
  kcdata.data.nfelem = std::get<2>(memsizes) + 0;
  kcdata.data.nielem = std::get<3>(memsizes) + 0;
  kcdata.data.nuelem = 2;
  kcdata.data.nbelem = 1;
}

template <typename T_IN, typename T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_device_(KernelConstantData<HGChefUncalibratedRecHitConstantData>& kcdata, cudautils::device::unique_ptr<double[]>& mem) {
  const std::vector<size_t> nelements = {kcdata.data.s_hgcHEF_fCPerMIP_, kcdata.data.s_hgcHEF_cce_, kcdata.data.s_hgcHEF_noise_fC_, kcdata.data.s_rcorr_, kcdata.data.s_weights_, kcdata.data.s_waferTypeL_};
  auto memsizes = get_memory_sizes_(nelements, 5, 0, 1);

  mem = cudautils::make_device_unique<double[]>(std::get<0>(memsizes), 0);

  kcdata.data.hgcHEF_fCPerMIP_ = mem.get();
  kcdata.data.hgcHEF_cce_      = kcdata.data.hgcHEF_fCPerMIP_ + nelements[0];
  kcdata.data.hgcHEF_noise_fC_ = kcdata.data.hgcHEF_cce_ + nelements[1];
  kcdata.data.rcorr_           = kcdata.data.hgcHEF_noise_fC_ + nelements[2];
  kcdata.data.weights_         = kcdata.data.rcorr_ + nelements[3];
  kcdata.data.waferTypeL_      = reinterpret_cast<int*>(kcdata.data.weights_ + nelements[4]);
  kcdata.data.nbytes = std::get<0>(memsizes);
  kcdata.data.ndelem = std::get<1>(memsizes) + 2;
  kcdata.data.nfelem = std::get<2>(memsizes) + 0;
  kcdata.data.nielem = std::get<3>(memsizes) + 0;
  kcdata.data.nuelem = 2;
  kcdata.data.nbelem = 1;
}


template <typename T_IN, typename T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_device_(KernelConstantData<HGChebUncalibratedRecHitConstantData>& kcdata, cudautils::device::unique_ptr<double[]>& mem) {
  const std::vector<size_t> nelements = {kcdata.data.s_hgcHEB_fCPerMIP_, kcdata.data.s_hgcHEB_cce_, kcdata.data.s_rcorr_, kcdata.data.s_weights_, kcdata.data.s_waferTypeL_};
  auto memsizes = get_memory_sizes_(nelements, 4, 0, 1);

  mem = cudautils::make_device_unique<double[]>(std::get<0>(memsizes), 0);

  kcdata.data.hgcHEB_fCPerMIP_ = mem.get();
  kcdata.data.hgcHEB_cce_      = kcdata.data.hgcHEB_fCPerMIP_ + nelements[0];
  kcdata.data.rcorr_           = kcdata.data.hgcHEB_cce_ + nelements[1];
  kcdata.data.weights_         = kcdata.data.rcorr_ + nelements[2];
  kcdata.data.waferTypeL_      = reinterpret_cast<int*>(kcdata.data.weights_ + nelements[3]);
  kcdata.data.nbytes = std::get<0>(memsizes);
  kcdata.data.ndelem = std::get<1>(memsizes) + 3;
  kcdata.data.nfelem = std::get<2>(memsizes) + 0;
  kcdata.data.nielem = std::get<3>(memsizes) + 0;
  kcdata.data.nuelem = 2;
  kcdata.data.nbelem = 1;
}

template <typename T_IN, typename T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_device_(HGCUncalibratedRecHitSoA*& soa1, HGCUncalibratedRecHitSoA*& soa2, HGCRecHitSoA*& soa3, cudautils::device::unique_ptr<float[]>& mem)
{
  std::vector<size_t> sizes = {6*sizeof(float), 3*sizeof(uint32_t),  //soa1
			       6*sizeof(float), 3*sizeof(uint32_t),  //soa2
			       2*sizeof(float), 3*sizeof(uint32_t)}; //soa3
  size_t size_tot = std::accumulate( sizes.begin(), sizes.end(), 0);
  mem = cudautils::make_device_unique<float[]>(stride_ * size_tot, 0);

  soa1->amplitude     = (float*)(mem.get());
  soa1->pedestal      = (float*)(soa1->amplitude    + stride_);
  soa1->jitter        = (float*)(soa1->pedestal     + stride_);
  soa1->chi2          = (float*)(soa1->jitter       + stride_);
  soa1->OOTamplitude  = (float*)(soa1->chi2         + stride_);
  soa1->OOTchi2       = (float*)(soa1->OOTamplitude + stride_);

  soa2->amplitude     = (float*)(soa1->OOTchi2      + stride_);
  soa2->pedestal      = (float*)(soa2->amplitude    + stride_);
  soa2->jitter        = (float*)(soa2->pedestal     + stride_);
  soa2->chi2          = (float*)(soa2->jitter       + stride_);
  soa2->OOTamplitude  = (float*)(soa2->chi2         + stride_);
  soa2->OOTchi2       = (float*)(soa2->OOTamplitude + stride_);

  soa3->energy        = (float*)(soa2->OOTchi2      + stride_);
  soa3->time          = (float*)(soa3->energy       + stride_);

  soa1->flags         = (uint32_t*)(soa3->time      + stride_);
  soa1->aux           = (uint32_t*)(soa1->flags     + stride_);
  soa1->id            = (uint32_t*)(soa1->aux       + stride_);

  soa2->flags         = (uint32_t*)(soa1->id        + stride_);
  soa2->aux           = (uint32_t*)(soa2->flags     + stride_);
  soa2->id            = (uint32_t*)(soa2->aux       + stride_);

  soa3->id            = (uint32_t*)(soa2->id        + stride_);
  soa3->flags         = (uint32_t*)(soa3->id        + stride_);
  soa3->flagBits      = (uint32_t*)(soa3->flags     + stride_);

  soa1->nbytes = std::accumulate(sizes.begin(), sizes.begin()+2, 0);
  soa2->nbytes = std::accumulate(sizes.begin()+2, sizes.begin()+4, 0);
  soa3->nbytes = std::accumulate(sizes.begin()+4, sizes.begin()+6, 0);
}

template <typename T_IN, typename T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_host_(KernelConstantData<HGCeeUncalibratedRecHitConstantData>& kcdata, cudautils::host::noncached::unique_ptr<double[]>& mem)
{
  const std::vector<size_t> nelements = {kcdata.data.s_hgcEE_fCPerMIP_, kcdata.data.s_hgcEE_cce_, kcdata.data.s_hgcEE_noise_fC_, kcdata.data.s_rcorr_, kcdata.data.s_weights_, kcdata.data.s_waferTypeL_};
  auto memsizes = get_memory_sizes_(nelements, 5, 0, 1);

  mem = cudautils::make_host_noncached_unique<double[]>(std::get<0>(memsizes), 0);

  kcdata.data.hgcEE_fCPerMIP_ = mem.get();
  kcdata.data.hgcEE_cce_      = kcdata.data.hgcEE_fCPerMIP_ + nelements[0];
  kcdata.data.hgcEE_noise_fC_ = kcdata.data.hgcEE_cce_ + nelements[1];
  kcdata.data.rcorr_          = kcdata.data.hgcEE_noise_fC_ + nelements[2];
  kcdata.data.weights_        = kcdata.data.rcorr_ + nelements[3];
  kcdata.data.waferTypeL_     = reinterpret_cast<int*>(kcdata.data.weights_ + nelements[4]);
  kcdata.data.nbytes = std::get<0>(memsizes);
  kcdata.data.ndelem = std::get<1>(memsizes) + 2;
  kcdata.data.nfelem = std::get<2>(memsizes) + 1;
  kcdata.data.nielem = std::get<3>(memsizes) + 0;
  kcdata.data.nuelem = 2;
  kcdata.data.nbelem = 1;
}

template <typename T_IN, typename T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_host_(KernelConstantData<HGChefUncalibratedRecHitConstantData>& kcdata, cudautils::host::noncached::unique_ptr<double[]>& mem)
{
  const std::vector<size_t> nelements = {kcdata.data.s_hgcHEF_fCPerMIP_, kcdata.data.s_hgcHEF_cce_, kcdata.data.s_hgcHEF_noise_fC_, kcdata.data.s_rcorr_, kcdata.data.s_weights_, kcdata.data.s_waferTypeL_};
  auto memsizes = get_memory_sizes_(nelements, 5, 0, 1);

  mem = cudautils::make_host_noncached_unique<double[]>(std::get<0>(memsizes), 0);

  kcdata.data.hgcHEF_fCPerMIP_ = mem.get();
  kcdata.data.hgcHEF_cce_      = kcdata.data.hgcHEF_fCPerMIP_ + nelements[0];
  kcdata.data.hgcHEF_noise_fC_ = kcdata.data.hgcHEF_cce_ + nelements[1];
  kcdata.data.rcorr_           = kcdata.data.hgcHEF_noise_fC_ + nelements[2];
  kcdata.data.weights_         = kcdata.data.rcorr_ + nelements[3];
  kcdata.data.waferTypeL_      = reinterpret_cast<int*>(kcdata.data.weights_ + nelements[4]);
  kcdata.data.nbytes = std::get<0>(memsizes);
  kcdata.data.ndelem = std::get<1>(memsizes) + 2;
  kcdata.data.nfelem = std::get<2>(memsizes) + 0;
  kcdata.data.nielem = std::get<3>(memsizes) + 0;
  kcdata.data.nuelem = 2;
  kcdata.data.nbelem = 1;
}

template <typename T_IN, typename T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_host_(KernelConstantData<HGChebUncalibratedRecHitConstantData>& kcdata, cudautils::host::noncached::unique_ptr<double[]>& mem)
{
  const std::vector<size_t> nelements = {kcdata.data.s_hgcHEB_fCPerMIP_, kcdata.data.s_hgcHEB_cce_, kcdata.data.s_rcorr_, kcdata.data.s_weights_, kcdata.data.s_waferTypeL_};
  auto memsizes = get_memory_sizes_(nelements, 4, 0, 1);

  mem = cudautils::make_host_noncached_unique<double[]>(std::get<0>(memsizes), 0);

  kcdata.data.hgcHEB_fCPerMIP_ = mem.get();
  kcdata.data.hgcHEB_cce_      = kcdata.data.hgcHEB_fCPerMIP_ + nelements[0];
  kcdata.data.rcorr_           = kcdata.data.hgcHEB_cce_ + nelements[1];
  kcdata.data.weights_         = kcdata.data.rcorr_ + nelements[2];
  kcdata.data.waferTypeL_      = reinterpret_cast<int*>(kcdata.data.weights_ + nelements[3]);
  kcdata.data.nbytes = std::get<0>(memsizes);
  kcdata.data.ndelem = std::get<1>(memsizes) + 3;
  kcdata.data.nfelem = std::get<2>(memsizes) + 0;
  kcdata.data.nielem = std::get<3>(memsizes) + 0;
  kcdata.data.nuelem = 2;
  kcdata.data.nbelem = 1;
}

template <typename T_IN, typename T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_host_(HGCUncalibratedRecHitSoA*& soa, cudautils::host::noncached::unique_ptr<float[]>& mem)
{
  size_t size1 = 6*sizeof(float);
  size_t size2 = 3*sizeof(uint32_t);
  mem = cudautils::make_host_noncached_unique<float[]>(stride_ * (size1+size2), 0);

  soa->amplitude     = (float*)(mem.get());
  soa->pedestal      = (float*)(soa->amplitude    + stride_);
  soa->jitter        = (float*)(soa->pedestal     + stride_);
  soa->chi2          = (float*)(soa->jitter       + stride_);
  soa->OOTamplitude  = (float*)(soa->chi2         + stride_);
  soa->OOTchi2       = (float*)(soa->OOTamplitude + stride_);
  soa->flags         = (uint32_t*)(soa->OOTchi2   + stride_);
  soa->aux           = (uint32_t*)(soa->flags     + stride_);
  soa->id            = (uint32_t*)(soa->aux       + stride_);
  soa->nbytes = size1 + size2;
}

template <typename T_IN, typename T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_host_(HGCRecHitSoA*& soa, cudautils::host::unique_ptr<float[]>& mem)
{
  size_t size1 = 2*sizeof(float);
  size_t size2 = 3*sizeof(uint32_t);
  mem = cudautils::make_host_unique<float[]>(stride_*(size1+size2), 0);

  soa->energy     = (float*)(mem.get());
  soa->time       = (float*)(soa->energy   + stride_);
  soa->id         = (uint32_t*)(soa->time  + stride_);
  soa->flags      = (uint32_t*)(soa->id    + stride_);
  soa->flagBits   = (uint32_t*)(soa->flags + stride_);
  soa->nbytes = size1 + size2;
}

template <typename T_IN, typename T_OUT>
template <class U>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::run(KernelConstantData<U>& kcdata)
{
  if (!std::is_same<T_IN, HGCUncalibratedRecHit>::value)
    throw cms::Exception("WrongTemplateType") << "The hgc_rechit_kernel_wrapper template does not support this type."; 
  set_geometry_constants_<U>(kcdata);

  KernelConstantData<U> d_kcdata = kcdata;
  //allocate pinned memory for constants on the host
  cudautils::host::noncached::unique_ptr<double[]> h_double_const;
  allocate_host_(kcdata, h_double_const);
  convert_constant_data_(kcdata);

  //allocate pinned memory for constants on the device
  cudautils::device::unique_ptr<double[]> d_double_const;
  allocate_device_(d_kcdata, d_double_const);

  //allocate memory for hits on the host
  old_soa_ = new HGCUncalibratedRecHitSoA();
  cudautils::host::noncached::unique_ptr<float[]> h_float_1;
  allocate_host_(old_soa_, h_float_1);
  convert_collection_data_to_soa_<HGCUncalibratedRecHit>();

  //allocate memory for hits on the device
  d_oldhits_ = new HGCUncalibratedRecHitSoA();
  d_newhits_ = new HGCUncalibratedRecHitSoA();
  d_newhits_final_ = new HGCRecHitSoA();
  cudautils::device::unique_ptr<float[]> d_float;
  allocate_device_(d_oldhits_, d_newhits_, d_newhits_final_, d_float);

  //allocate memory for hits on the host
  h_newhits_ = new HGCRecHitSoA();
  cudautils::host::unique_ptr<float[]> h_float_2;
  allocate_host_(h_newhits_, h_float_2);

  KernelModifiableData<HGCUncalibratedRecHitSoA, HGCRecHitSoA> kmdata(nhits_, old_soa_, d_oldhits_, d_newhits_, d_newhits_final_, h_newhits_);
  KernelManagerHGCalRecHit kernel_manager(kmdata, det_);

  kernel_manager.run_kernels(kcdata, d_kcdata);
  new_soa_ = kernel_manager.get_output();
  convert_soa_data_to_collection_<HGCRecHit>();
}

template <class T_IN, class T_OUT>
edm::SortedCollection<T_OUT> HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::get_output_collection() 
{
  return out_data_;
}

template <class T_IN, class T_OUT>
template <class U_IN>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::convert_collection_data_to_soa_() 
{
  throw cms::Exception("NoSpecialization") << "convert_collection_data_to_soa_()";  
}

template <class T_IN, class T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::convert_constant_data_(KernelConstantData<HGCeeUncalibratedRecHitConstantData>& kcdata)
{
  for(unsigned int i=0; i<kcdata.data.s_hgcEE_fCPerMIP_; ++i)
    kcdata.data.hgcEE_fCPerMIP_[i] = kcdata.vdata.fCPerMIP[i];
  for(unsigned int i=0; i<kcdata.data.s_hgcEE_cce_; ++i)
    kcdata.data.hgcEE_cce_[i] = kcdata.vdata.cce[i];
  for(unsigned int i=0; i<kcdata.data.s_hgcEE_noise_fC_; ++i)
    kcdata.data.hgcEE_noise_fC_[i] = kcdata.vdata.noise_fC[i];
  for(unsigned int i=0; i<kcdata.data.s_rcorr_; ++i)
    kcdata.data.rcorr_[i] = kcdata.vdata.rcorr[i];
  for(unsigned int i=0; i<kcdata.data.s_weights_; ++i)
    kcdata.data.weights_[i] = kcdata.vdata.weights[i];
  for(unsigned int i=0; i<kcdata.data.s_waferTypeL_; ++i)
    kcdata.data.waferTypeL_[i] = kcdata.vdata.waferTypeL[i];
}

template <class T_IN, class T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::convert_constant_data_(KernelConstantData<HGChefUncalibratedRecHitConstantData>& kcdata)
{
  for(unsigned int i=0; i<kcdata.data.s_hgcHEF_fCPerMIP_; ++i)
    kcdata.data.hgcHEF_fCPerMIP_[i] = kcdata.vdata.fCPerMIP[i];
  for(unsigned int i=0; i<kcdata.data.s_hgcHEF_cce_; ++i)
    kcdata.data.hgcHEF_cce_[i] = kcdata.vdata.cce[i];
  for(unsigned int i=0; i<kcdata.data.s_hgcHEF_noise_fC_; ++i)
    kcdata.data.hgcHEF_noise_fC_[i] = kcdata.vdata.noise_fC[i];
  for(unsigned int i=0; i<kcdata.data.s_rcorr_; ++i)
    kcdata.data.rcorr_[i] = kcdata.vdata.rcorr[i];
  for(unsigned int i=0; i<kcdata.data.s_weights_; ++i)
    kcdata.data.weights_[i] = kcdata.vdata.weights[i];
  for(unsigned int i=0; i<kcdata.data.s_waferTypeL_; ++i)
    kcdata.data.waferTypeL_[i] = kcdata.vdata.waferTypeL[i];
}

template <class T_IN, class T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::convert_constant_data_(KernelConstantData<HGChebUncalibratedRecHitConstantData>& kcdata)
{
  for(unsigned int i=0; i<kcdata.data.s_hgcHEB_fCPerMIP_; ++i)
    kcdata.data.hgcHEB_fCPerMIP_[i] = kcdata.vdata.fCPerMIP[i];
  for(unsigned int i=0; i<kcdata.data.s_hgcHEB_cce_; ++i)
    kcdata.data.hgcHEB_cce_[i] = kcdata.vdata.cce[i];
  for(unsigned int i=0; i<kcdata.data.s_rcorr_; ++i)
    kcdata.data.rcorr_[i] = kcdata.vdata.rcorr[i];
  for(unsigned int i=0; i<kcdata.data.s_weights_; ++i)
    kcdata.data.weights_[i] = kcdata.vdata.weights[i];
  for(unsigned int i=0; i<kcdata.data.s_waferTypeL_; ++i)
    kcdata.data.waferTypeL_[i] = kcdata.vdata.waferTypeL[i];
}

template <>
template <> 
void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_collection_data_to_soa_<HGCUncalibratedRecHit>()
{
  for(unsigned int i=0; i<nhits_; ++i)
    {
      old_soa_->amplitude[i] = hits_[i].amplitude();
      old_soa_->pedestal[i] = hits_[i].pedestal();
      old_soa_->jitter[i] = hits_[i].jitter();
      old_soa_->chi2[i] = hits_[i].chi2();
      old_soa_->OOTamplitude[i] = hits_[i].outOfTimeEnergy();
      old_soa_->OOTchi2[i] = hits_[i].outOfTimeChi2();
      old_soa_->flags[i] = hits_[i].flags();
      old_soa_->aux[i] = 0;
      old_soa_->id[i] = hits_[i].id().rawId();
    }
}

template <class T_IN, class T_OUT>
template <class U_IN>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::convert_soa_data_to_collection_() 
{
  throw cms::Exception("NoSpecialization") << "convert_soa_data_to_collection_";
}

template <>
template <> 
void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_soa_data_to_collection_<HGCRecHit>()
{
  for(uint i=0; i<nhits_; ++i)
    {
      DetId id_converted( new_soa_->id[i] );
      out_data_[i] = HGCRecHit(id_converted, new_soa_->energy[i], new_soa_->time[i], new_soa_->flags[i], new_soa_->flagBits[i]);
    }
}

template class HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>;
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::set_geometry_constants_<HGCeeUncalibratedRecHitConstantData>(KernelConstantData<HGCeeUncalibratedRecHitConstantData>&);
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::set_geometry_constants_<HGChefUncalibratedRecHitConstantData>(KernelConstantData<HGChefUncalibratedRecHitConstantData>&);
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::set_geometry_constants_<HGChebUncalibratedRecHitConstantData>(KernelConstantData<HGChebUncalibratedRecHitConstantData>&);
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::run<HGCeeUncalibratedRecHitConstantData>(KernelConstantData<HGCeeUncalibratedRecHitConstantData>&);
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::run<HGChefUncalibratedRecHitConstantData>(KernelConstantData<HGChefUncalibratedRecHitConstantData>&);
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::run<HGChebUncalibratedRecHitConstantData>(KernelConstantData<HGChebUncalibratedRecHitConstantData>&);
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_collection_data_to_soa_<HGCUncalibratedRecHit>();
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_soa_data_to_collection_<HGCRecHit>();
