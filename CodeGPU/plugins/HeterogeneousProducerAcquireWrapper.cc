#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>
#include <type_traits>
#include <cuda_runtime.h>

#include "HeterogeneousProducerAcquireWrapper.h"
#include "KernelManager.h"
#include "Types.h"

/*
I am currently mixing GPU and non-GPU types 
Eventually all types will be given by T_IN and T_OUT
*/
template <class T_IN, class T_OUT>
HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::HeterogeneousProducerAcquireWrapper(const edm::SortedCollection<T_IN>& hits, const edm::EventSetup& setup)
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
HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::~HeterogeneousProducerAcquireWrapper()
{
}

template <class T_IN, class T_OUT>
void HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::set_geometry_(const edm::EventSetup& setup)
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
void HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::allocate_device_(U& soa, cudautils::device::unique_ptr<float[]>& mem)
{
  throw cms::Exception("NoSpecialization") << "allocate_device_()";  
}

//specialization for HGCUncalibratedRecHitSoA
template <>
template <> 
void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_device_<HGCUncalibratedRecHitSoA>(HGCUncalibratedRecHitSoA& soa, cudautils::device::unique_ptr<float[]>& mem)
{
  const size_t size1 = 6*stride_*sizeof(float);
  const size_t size2 = 3*stride_*sizeof(uint32_t);
  mem = cudautils::make_device_unique<float[]>(size1 + size2, 0);
  soa.amplitude     = mem.get();
  soa.pedestal      = mem.get() +   stride_;
  soa.jitter        = mem.get() + 2*stride_;
  soa.chi2          = mem.get() + 3*stride_;
  soa.OOTamplitude  = mem.get()    + 4*stride_;
  soa.OOTchi2       = mem.get()    + 5*stride_;
  soa.flags         = (uint32_t*)mem.get() + 6*stride_;
  soa.aux           = (uint32_t*)mem.get() + 7*stride_;
  soa.id            = (uint32_t*)mem.get() + 8*stride_;
  soa.nbytes = size1 + size2;
}

//specialization for HGCRecHitSoA
template <>
template <> 
void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_device_<HGCRecHitSoA>(HGCRecHitSoA& soa, cudautils::device::unique_ptr<float[]>& mem)
{
  const size_t size1 = 2*stride_*sizeof(float);
  const size_t size2 = 3*stride_*sizeof(uint32_t);
  mem = cudautils::make_device_unique<float[]>(size1 + size2, 0); 
  soa.energy     = mem.get();
  soa.time       = mem.get() +   stride_;
  soa.id         = (uint32_t*)mem.get() + 2*stride_;
  soa.flags      = (uint32_t*)mem.get() + 3*stride_;
  soa.flagBits   = (uint32_t*)mem.get() + 4*stride_;
  soa.nbytes = size1 + size2;
}


template <typename T_IN, typename T_OUT>
template <class U>
void HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::allocate_host_(U& soa, cudautils::host::unique_ptr<float[]>& mem)
{
  throw cms::Exception("NoSpecialization") << "allocate_host_() cached";  
}

//overload for non-cached (pinned/non-pageable) memory
template <typename T_IN, typename T_OUT>
template <class U>
void HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::allocate_host_(U& soa, cudautils::host::noncached::unique_ptr<float[]>& mem)
{
  throw cms::Exception("NoSpecialization") << "allocate_host_() non-cached";  
}

//specialization for HGCUncalibratedRecHitSoA
template <>
template <> 
void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_host_<HGCUncalibratedRecHitSoA>(HGCUncalibratedRecHitSoA& soa, cudautils::host::noncached::unique_ptr<float[]>& mem)
{
  const unsigned int size1 = 6*stride_*sizeof(float);
  const unsigned int size2 = 3*stride_*sizeof(uint32_t);
  mem = cudautils::make_host_noncached_unique<float[]>(size1 + size2, 0);
  soa.amplitude     = mem.get();
  soa.pedestal      = mem.get()    +   stride_;
  soa.jitter        = mem.get()    + 2*stride_;
  soa.chi2          = mem.get()    + 3*stride_;
  soa.OOTamplitude  = mem.get()    + 4*stride_;
  soa.OOTchi2       = mem.get()    + 5*stride_;
  soa.flags         = (uint32_t*)mem.get()   + 6*stride_;
  soa.aux           = (uint32_t*)mem.get()   + 7*stride_;
  soa.id            = (uint32_t*)mem.get()   + 8*stride_;
  soa.nbytes = size1 + size2;
}

//specialization for HGCRecHitSoA
template <>
template <> 
void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_host_<HGCRecHitSoA>(HGCRecHitSoA& soa, cudautils::host::unique_ptr<float[]>& mem)
{
  const unsigned int size1 = 2*stride_*sizeof(float);
  const unsigned int size2 = 3*stride_*sizeof(uint32_t);
  mem = cudautils::make_host_unique<float[]>(size1 + size2, 0);
  soa.energy     = mem.get();
  soa.time       = mem.get()  +  stride_;
  soa.id         = (uint32_t*)mem.get() + 2*stride_;
  soa.flags      = (uint32_t*)mem.get() + 3*stride_;
  soa.flagBits   = (uint32_t*)mem.get() + 4*stride_;
  soa.nbytes = size1 + size2;
}

template <typename T_IN, typename T_OUT>
void HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::run()
{
  if (!std::is_same<T_IN, HGCUncalibratedRecHit>::value)
    throw cms::Exception("WrongTemplateType") << "The hgc_rechit_kernel_wrapper template does not support this type.";
  
  HGCUncalibratedRecHitSoA old_soa_obj;
  cudautils::host::noncached::unique_ptr<float[]> h_float_1;
  allocate_host_<HGCUncalibratedRecHitSoA>(old_soa_obj, h_float_1);
  HGCUncalibratedRecHitSoA *old_soa = &old_soa_obj;

  convert_collection_data_to_soa_<HGCUncalibratedRecHit>();

  HGCUncalibratedRecHitSoA d_soa_2;
  cudautils::device::unique_ptr<float[]> d_float_2;
  allocate_device_<HGCUncalibratedRecHitSoA>(d_soa_2, d_float_2);
  HGCUncalibratedRecHitSoA *d_oldhits = &d_soa_2;

  HGCUncalibratedRecHitSoA d_soa_3;
  cudautils::device::unique_ptr<float[]> d_float_3;
  allocate_device_<HGCUncalibratedRecHitSoA>(d_soa_3, d_float_3);
  HGCUncalibratedRecHitSoA *d_newhits = &d_soa_3;

  HGCRecHitSoA d_soa_4;
  cudautils::device::unique_ptr<float[]> d_float_4;
  allocate_device_<HGCRecHitSoA>(d_soa_4, d_float_4);
  HGCRecHitSoA *d_newhits_final = &d_soa_4;

  HGCRecHitSoA h_soa_2;
  cudautils::host::unique_ptr<float[]> h_float_2;
  allocate_host_<HGCRecHitSoA>(h_soa_2, h_float_2);
  HGCRecHitSoA *h_newhits = &h_soa_2;

  assert(old_soa->nbytes == d_oldhits->nbytes);
  assert(old_soa->nbytes == d_newhits->nbytes);
  assert(d_newhits_final->nbytes == h_newhits->nbytes);
  KernelManagerData<HGCUncalibratedRecHitSoA, HGCRecHitSoA> kmdata(nhits_, 
				  old_soa, d_oldhits, d_newhits, 
                                  d_newhits_final, h_newhits);

  DetId::Detector dtype = DetId::Detector::HGCalEE;
  KernelManagerHGCalRecHit kernel_manager(kmdata, dtype);
  std::cout << "before run_kernels()" << std::endl;
  kernel_manager.run_kernels();
  std::cout << "after run_kernels()" << std::endl;
  new_soa_ = kernel_manager.get_output();

  convert_soa_data_to_collection_<HGCRecHit>();
}

template <class T_IN, class T_OUT>
edm::SortedCollection<T_OUT> HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::get_output_collection() 
{
  return out_data_;
}

template <class T_IN, class T_OUT>
template <class U_IN>
void HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::convert_collection_data_to_soa_() 
{
  throw cms::Exception("NoSpecialization") << "convert_collection_data_to_soa_()";  
}

template <>
template <> 
void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_collection_data_to_soa_<HGCUncalibratedRecHit>()
{
  for(unsigned int i=0; i<nhits_; ++i)
    {
      std::cout << hits_[i].amplitude() << std::endl;
      /*
      old_soa_->amplitude[i] = hits_[i].amplitude();
      old_soa_->pedestal[i] = hits_[i].pedestal();
      old_soa_->jitter[i] = hits_[i].jitter();
      old_soa_->chi2[i] = hits_[i].chi2();
      old_soa_->OOTamplitude[i] = hits_[i].outOfTimeEnergy();
      old_soa_->OOTchi2[i] = hits_[i].outOfTimeChi2();
      old_soa_->flags[i] = hits_[i].flags();
      old_soa_->aux[i] = 0;
      old_soa_->id[i] = hits_[i].id().rawId();
*/
    }
}

template <class T_IN, class T_OUT>
template <class U_IN>
void HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::convert_soa_data_to_collection_() 
{
  throw cms::Exception("NoSpecialization") << "convert_soa_data_to_collection_";
}

template <>
template <> 
void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_soa_data_to_collection_<HGCRecHit>()
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

template class HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>;
template void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_collection_data_to_soa_<HGCUncalibratedRecHit>();
template void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_soa_data_to_collection_<HGCRecHit>();

template void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_host_<HGCUncalibratedRecHitSoA>(HGCUncalibratedRecHitSoA&, cudautils::host::noncached::unique_ptr<float[]>&);
template void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_host_<HGCRecHitSoA>(HGCRecHitSoA&, cudautils::host::unique_ptr<float[]>&);
template void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_device_<HGCUncalibratedRecHitSoA>(HGCUncalibratedRecHitSoA&, cudautils::device::unique_ptr<float[]>&);
template void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::allocate_device_<HGCRecHitSoA>(HGCRecHitSoA&, cudautils::device::unique_ptr<float[]>&);
