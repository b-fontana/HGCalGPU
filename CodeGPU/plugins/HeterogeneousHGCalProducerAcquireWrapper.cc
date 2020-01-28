#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>
#include <type_traits>
#include <cuda_runtime.h>

#include "HeterogeneousHGCalProducerAcquireWrapper.h"
#include "KernelManager.h"
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
  if( (det_ != DetId::HGCalEE) and (det_ != DetId::HGCalHSi) and (det_ != DetId::HGCalHSc))
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

template <typename T_IN, typename T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_device_(HGCUncalibratedRecHitSoA*& soa, cudautils::device::unique_ptr<float[]>& mem)
{
  const size_t size1 = 6*stride_*sizeof(float);
  const size_t size2 = 3*stride_*sizeof(uint32_t);
  mem = cudautils::make_device_unique<float[]>(size1 + size2, 0);

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
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_device_(HGCRecHitSoA*& soa, cudautils::device::unique_ptr<float[]>& mem)
{
  const size_t size1 = 2*stride_*sizeof(float);
  const size_t size2 = 3*stride_*sizeof(uint32_t);
  mem = cudautils::make_device_unique<float[]>(size1 + size2, 0); 

  soa->energy     = (float*)(mem.get());
  soa->time       = (float*)(soa->energy   + stride_);
  soa->id         = (uint32_t*)(soa->time  + stride_);
  soa->flags      = (uint32_t*)(soa->id    + stride_);
  soa->flagBits   = (uint32_t*)(soa->flags + stride_);
  soa->nbytes = size1 + size2;
}

template <typename T_IN, typename T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_host_(HGCUncalibratedRecHitSoA*& soa, cudautils::host::noncached::unique_ptr<float[]>& mem)
{
  const size_t size1 = 6*stride_*sizeof(float);
  const size_t size2 = 3*stride_*sizeof(uint32_t);
  mem = cudautils::make_host_noncached_unique<float[]>(size1 + size2, 0);

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
  const size_t size1 = 2*stride_*sizeof(float);
  const size_t size2 = 3*stride_*sizeof(uint32_t);
  mem = cudautils::make_host_unique<float[]>(size1 + size2, 0);

  soa->energy     = (float*)(mem.get());
  soa->time       = (float*)(soa->energy   + stride_);
  soa->id         = (uint32_t*)(soa->time  + stride_);
  soa->flags      = (uint32_t*)(soa->id    + stride_);
  soa->flagBits   = (uint32_t*)(soa->flags + stride_);
  soa->nbytes = size1 + size2;
}

template <typename T_IN, typename T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::run()
{
  if (!std::is_same<T_IN, HGCUncalibratedRecHit>::value)
    throw cms::Exception("WrongTemplateType") << "The hgc_rechit_kernel_wrapper template does not support this type.";

  old_soa_ = new HGCUncalibratedRecHitSoA();
  cudautils::host::noncached::unique_ptr<float[]> h_float_1;
  allocate_host_(old_soa_, h_float_1);
  convert_collection_data_to_soa_<HGCUncalibratedRecHit>();
  d_oldhits_ = new HGCUncalibratedRecHitSoA();
  cudautils::device::unique_ptr<float[]> d_float_2;
  std::cout << "check0 2"  << std::endl;
  allocate_device_(d_oldhits_, d_float_2);
  std::cout << "check1" << std::endl;
  d_newhits_ = new HGCUncalibratedRecHitSoA();
  cudautils::device::unique_ptr<float[]> d_float_3;
  allocate_device_(d_newhits_, d_float_3);
  std::cout << "check2" << std::endl;
  d_newhits_final_ = new HGCRecHitSoA();
  cudautils::device::unique_ptr<float[]> d_float_4;
  allocate_device_(d_newhits_final_, d_float_4);
  std::cout << "check3" << std::endl;
  h_newhits_ = new HGCRecHitSoA();
  cudautils::host::unique_ptr<float[]> h_float_2;
  allocate_host_(h_newhits_, h_float_2);
  std::cout << "check4" << std::endl;
  /*
  assert(old_soa_->nbytes == d_oldhits_->nbytes);
  assert(old_soa_->nbytes == d_newhits_->nbytes);
  assert(d_newhits_final_->nbytes == h_newhits_->nbytes);
  */
  std::cout << "check5" << std::endl;
  KernelManagerData<HGCUncalibratedRecHitSoA, HGCRecHitSoA> kmdata(nhits_, old_soa_, d_oldhits_, d_newhits_, d_newhits_final_, h_newhits_);
  std::cout << "check6" << std::endl;
  DetId::Detector dtype = DetId::Detector::HGCalEE;
  KernelManagerHGCalRecHit kernel_manager(kmdata, dtype);
  std::cout << "before run_kernels()" << std::endl;
  kernel_manager.run_kernels();
  std::cout << "after run_kernels()" << std::endl;
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
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_collection_data_to_soa_<HGCUncalibratedRecHit>();
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_soa_data_to_collection_<HGCRecHit>();
