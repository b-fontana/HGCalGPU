#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>
#include <type_traits>
#include <cuda_runtime.h>

#include "HeterogeneousHGCalProducerAcquireWrapper.h"
#include "KernelManager.h"
#include "Types.h"

/*
I am currently mixing GPU and non-GPU types 
Eventually all types will be given by T_IN and T_OUT
*/
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
template <class U>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_device_(U* soa, cudautils::device::unique_ptr<U[]>& mem)
{
  throw cms::Exception("NoSpecialization") << "allocate_device_()";  
}

//specialization for HGCUncalibratedRecHitSoA
template <>
template <> 
void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_device_<HGCUncalibratedRecHitSoA>(HGCUncalibratedRecHitSoA* soa, cudautils::device::unique_ptr<HGCUncalibratedRecHitSoA[]>& mem)
{
  const size_t size1 = sizeof(HGCUncalibratedRecHitSoA);
  const size_t size2 = 6*stride_*sizeof(float);
  const size_t size3 = 3*stride_*sizeof(uint32_t);
  mem = cudautils::make_device_unique<HGCUncalibratedRecHitSoA[]>(size1 + size2 + size3, 0);
  soa                = mem.get();
  soa->amplitude     = (float*)(soa + 1);
  soa->pedestal      = (float*)(soa->amplitude    + stride_);
  soa->jitter        = (float*)(soa->pedestal     + stride_);
  soa->chi2          = (float*)(soa->jitter       + stride_);
  soa->OOTamplitude  = (float*)(soa->chi2         + stride_);
  soa->OOTchi2       = (float*)(soa->OOTamplitude + stride_);
  soa->flags         = (uint32_t*)(soa->OOTchi2   + stride_);
  soa->aux           = (uint32_t*)(soa->flags     + stride_);
  soa->id            = (uint32_t*)(soa->aux       + stride_);
  soa->nbytes = size1 + size2 + size3;
}

//specialization for HGCRecHitSoA
template <>
template <> 
void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_device_<HGCRecHitSoA>(HGCRecHitSoA* soa, cudautils::device::unique_ptr<HGCRecHitSoA[]>& mem)
{
  const size_t size1 = sizeof(HGCRecHitSoA);
  const size_t size2 = 2*stride_*sizeof(float);
  const size_t size3 = 3*stride_*sizeof(uint32_t);
  mem = cudautils::make_device_unique<HGCRecHitSoA[]>(size1 + size2 + size3, 0); 
  soa            = mem.get();
  soa->energy     = (float*)mem.get() + stride_;
  soa->time       = (float*)mem.get() + 2*stride_;
  soa->id         = (uint32_t*)mem.get() + 3*stride_;
  soa->flags      = (uint32_t*)mem.get() + 4*stride_;
  soa->flagBits   = (uint32_t*)mem.get() + 5*stride_;
  soa->nbytes = size1 + size2 + size3;
}


template <typename T_IN, typename T_OUT>
template <class U>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_host_(U* soa, cudautils::host::unique_ptr<U[]>& mem)
{
  throw cms::Exception("NoSpecialization") << "allocate_host_() cached";  
}

//overload for non-cached (pinned/non-pageable) memory
template <typename T_IN, typename T_OUT>
template <class U>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::allocate_host_(U* soa, cudautils::host::noncached::unique_ptr<U[]>& mem)
{
  throw cms::Exception("NoSpecialization") << "allocate_host_() non-cached";  
}

//specialization for HGCUncalibratedRecHitSoA
template <>
template <> 
void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_host_<HGCUncalibratedRecHitSoA>(HGCUncalibratedRecHitSoA* soa, cudautils::host::noncached::unique_ptr<HGCUncalibratedRecHitSoA[]>& mem)
{
  const size_t size1 = sizeof(HGCUncalibratedRecHitSoA);
  const size_t size2 = 6*stride_*sizeof(float);
  const size_t size3 = 3*stride_*sizeof(uint32_t);
  mem = cudautils::make_host_noncached_unique<HGCUncalibratedRecHitSoA[]>(size1 + size2 + size3, 0);
  soa                = mem.get();
  soa->amplitude     = (float*)(soa + 1);
  soa->pedestal      = (float*)(soa->amplitude    + stride_);
  soa->jitter        = (float*)(soa->pedestal     + stride_);
  soa->chi2          = (float*)(soa->jitter       + stride_);
  soa->OOTamplitude  = (float*)(soa->chi2         + stride_);
  soa->OOTchi2       = (float*)(soa->OOTamplitude + stride_);
  soa->flags         = (uint32_t*)(soa->OOTchi2   + stride_);
  soa->aux           = (uint32_t*)(soa->flags     + stride_);
  soa->id            = (uint32_t*)(soa->aux       + stride_);
  soa->nbytes = size1 + size2 + size3;
}

//specialization for HGCRecHitSoA
template <>
template <> 
void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_host_<HGCRecHitSoA>(HGCRecHitSoA* soa, cudautils::host::unique_ptr<HGCRecHitSoA[]>& mem)
{
  const size_t size1 = sizeof(HGCRecHitSoA);
  const size_t size2 = 2*stride_*sizeof(float);
  const size_t size3 = 3*stride_*sizeof(uint32_t);
  mem = cudautils::make_host_unique<HGCRecHitSoA[]>(size1 + size2 + size3, 0);
  soa            = mem.get();
  soa->energy     = (float*)(soa + 1);
  soa->time       = (float*)(mem.get() + 2*stride_;
  soa->id         = (uint32_t*)mem.get() + 3*stride_;
  soa->flags      = (uint32_t*)mem.get() + 4*stride_;
  soa->flagBits   = (uint32_t*)mem.get() + 5*stride_;
  soa->nbytes = size1 + size2 + size3;
}

template <typename T_IN, typename T_OUT>
void HeterogeneousHGCalProducerAcquireWrapper<T_IN, T_OUT>::run()
{
  if (!std::is_same<T_IN, HGCUncalibratedRecHit>::value)
    throw cms::Exception("WrongTemplateType") << "The hgc_rechit_kernel_wrapper template does not support this type.";

  std::cout << "check" << std::endl;
  cudautils::host::noncached::unique_ptr<HGCUncalibratedRecHitSoA[]> h_float_1;
  allocate_host_<HGCUncalibratedRecHitSoA>(old_soa_, h_float_1);

  convert_collection_data_to_soa_<HGCUncalibratedRecHit>();
  std::cout << "check0" << std::endl;
  HGCUncalibratedRecHitSoA *d_oldhits = nullptr;
  cudautils::device::unique_ptr<HGCUncalibratedRecHitSoA[]> d_float_2;
  allocate_device_<HGCUncalibratedRecHitSoA>(d_oldhits, d_float_2);
  std::cout << "check1" << std::endl;
  HGCUncalibratedRecHitSoA *d_newhits = nullptr;
  cudautils::device::unique_ptr<HGCUncalibratedRecHitSoA[]> d_float_3;
  allocate_device_<HGCUncalibratedRecHitSoA>(d_newhits, d_float_3);
  std::cout << "check2" << std::endl;
  HGCRecHitSoA *d_newhits_final = nullptr;
  cudautils::device::unique_ptr<HGCRecHitSoA[]> d_float_4;
  allocate_device_<HGCRecHitSoA>(d_newhits_final, d_float_4);
  std::cout << "check3" << std::endl;
  HGCRecHitSoA *h_newhits = nullptr;
  cudautils::host::unique_ptr<HGCRecHitSoA[]> h_float_2;
  allocate_host_<HGCRecHitSoA>(h_newhits, h_float_2);
  std::cout << "check4" << std::endl;
  assert(old_soa_->nbytes == d_oldhits->nbytes);
  assert(old_soa_->nbytes == d_newhits->nbytes);
  assert(d_newhits_final->nbytes == h_newhits->nbytes);
  std::cout << "check5" << std::endl;
  KernelManagerData<HGCUncalibratedRecHitSoA, HGCRecHitSoA> kmdata(nhits_, old_soa_, d_oldhits, d_newhits, d_newhits_final, h_newhits);
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
  for(uint i=0; i<10/*nhits_*/; ++i)
    {
      DetId id_converted( new_soa_->id[i] );
      out_data_[i] = HGCRecHit(id_converted,
			       new_soa_->energy[i], 
			       new_soa_->time[i], 
			       new_soa_->flags[i], 
			       new_soa_->flagBits[i]);
    }
}

template class HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>;
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_collection_data_to_soa_<HGCUncalibratedRecHit>();
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_soa_data_to_collection_<HGCRecHit>();

template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_host_<HGCUncalibratedRecHitSoA>(HGCUncalibratedRecHitSoA*, cudautils::host::noncached::unique_ptr<HGCUncalibratedRecHitSoA[]>&);
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_host_<HGCRecHitSoA>(HGCRecHitSoA*, cudautils::host::unique_ptr<HGCRecHitSoA[]>&);
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_device_<HGCUncalibratedRecHitSoA>(HGCUncalibratedRecHitSoA*, cudautils::device::unique_ptr<HGCUncalibratedRecHitSoA[]>&);
template void HeterogeneousHGCalProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_device_<HGCRecHitSoA>(HGCRecHitSoA*, cudautils::device::unique_ptr<HGCRecHitSoA[]>&);
