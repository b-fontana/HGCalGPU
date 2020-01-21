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
  size_ = hits.size();
  if (size_ == 0)
    throw cms::Exception("EmptyCollection") << "The passed collection is empty."; 
 
  for(unsigned int i=0; i<size_-1; ++i)
    assert(hits[i].id().det() == hits[i+1].id().det());
  
  det_ = hits[0].id().det(); 
  if( (det_ != DetId::HGCalEE) and (det_ != DetId::HGCalHSi) and (det_ != DetId::HGCalHSc))
    throw cms::Exception("WrongDetectorType") << "The specified detector is wrong.";

  hits_ = hits;
  old_soa_ = new HGCUncalibratedRecHitSoA;
  new_soa_ = new HGCRecHitSoA;
  out_data_.reserve(size_);
  tools_.reset(new hgcal::RecHitTools());
  this->set_geometry_(setup);
}

template <class T_IN, class T_OUT>
HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::~HeterogeneousProducerAcquireWrapper()
{
  delete old_soa_->amplitude;     old_soa_->amplitude    = nullptr;
  delete old_soa_->pedestal;      old_soa_->pedestal     = nullptr;
  delete old_soa_->jitter;        old_soa_->jitter       = nullptr;
  delete old_soa_->chi2;          old_soa_->chi2         = nullptr;
  delete old_soa_->OOTamplitude;  old_soa_->OOTamplitude = nullptr;
  delete old_soa_->OOTchi2;       old_soa_->OOTchi2      = nullptr;
  delete old_soa_->flags;         old_soa_->flags        = nullptr;
  delete old_soa_->aux;           old_soa_->aux          = nullptr;
  delete old_soa_->id;            old_soa_->id           = nullptr;

  delete old_soa_;                old_soa_ = nullptr;
  delete new_soa_;                new_soa_ = nullptr;
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
void HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::run()
{
  if (!std::is_same<T_IN, HGCUncalibratedRecHit>::value)
    throw cms::Exception("WrongTemplateType") << "The hgc_rechit_kernel_wrapper template does not support this type.";
  
  convert_collection_data_to_soa_<HGCUncalibratedRecHit>();

  auto cpu_input = cudautils::make_host_noncached_unique< HGCUncalibratedRecHitSoA[] >(size_, cudaHostAllocWriteCombined);
  std::memcpy(cpu_input.get(), old_soa_, sizeof(HGCUncalibratedRecHitSoA) * size_);
  auto d_oldhits = cudautils::make_device_unique<HGCUncalibratedRecHitSoA[]>(size_, 0);
  auto d_newhits = cudautils::make_device_unique<HGCUncalibratedRecHitSoA[]>(size_, 0);
  auto d_newhits_final = cudautils::make_device_unique<HGCRecHitSoA[]>(size_, 0);
  auto h_newhits = cudautils::make_host_unique<HGCRecHitSoA[]>(size_, 0);
  KernelManagerData<HGCUncalibratedRecHitSoA, HGCRecHitSoA> kmdata(cpu_input.get(), d_oldhits.get(), d_newhits.get(), d_newhits_final.get(), h_newhits.get());
  DetId::Detector dtype = DetId::Detector::HGCalEE;
  KernelManagerHGCalRecHit kernel_manager(kmdata, dtype);
  kernel_manager.run_kernels();
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
}

template <>
template <> 
void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_collection_data_to_soa_<HGCUncalibratedRecHit>()
{
  old_soa_->size = size_;

  old_soa_->amplitude    = new float    [ old_soa_->size ];
  old_soa_->pedestal     = new float    [ old_soa_->size ];
  old_soa_->jitter       = new float    [ old_soa_->size ];
  old_soa_->chi2         = new float    [ old_soa_->size ];
  old_soa_->OOTamplitude = new float    [ old_soa_->size ];
  old_soa_->OOTchi2      = new float    [ old_soa_->size ];
  old_soa_->flags        = new uint32_t [ old_soa_->size ];
  old_soa_->aux          = new uint32_t [ old_soa_->size ];
  old_soa_->id           = new uint32_t [ old_soa_->size ];

  for(uint i=0; i<old_soa_->size; ++i)
    {
      old_soa_->amplitude[i]    = hits_[i].amplitude();
      old_soa_->pedestal[i]     = hits_[i].pedestal();
      old_soa_->jitter[i]       = hits_[i].jitter();
      old_soa_->chi2[i]         = hits_[i].chi2();
      old_soa_->OOTamplitude[i] = hits_[i].outOfTimeEnergy();
      old_soa_->OOTchi2[i]      = hits_[i].outOfTimeChi2();
      old_soa_->flags[i]        = hits_[i].flags();
      old_soa_->aux[i]          = hits_[i].jitterErrorBits();
      old_soa_->id[i]           = hits_[i].id().rawId();
    }
}

template <class T_IN, class T_OUT>
template <class U_IN>
void HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::convert_soa_data_to_collection_() 
{
}

template <>
template <> 
void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_soa_data_to_collection_<HGCRecHit>()
{
  for(uint i=0; i<10/*size_*/; ++i)
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
