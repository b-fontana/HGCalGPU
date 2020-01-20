#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>
#include <type_traits>
#include <cuda_runtime.h>

#include "HeterogeneousProducerAcquireWrapper.h"
#include "KernelManager.h"
#include "HGCalRecHitKernelImpl.cuh"

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
  tools_.reset(new hgcal::RecHitTools());
  this->set_geometry_(setup);
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

  //eventually this data will be produced in a 'convert_collection_data_to_soa_()'
  std::vector<HGCUncalibratedRecHit_GPU> aux_in(size_);
  for(uint i=0; i<size_; ++i)
    {
      HGCUncalibratedRecHit_GPU hit;
      aux_in[i] = hit;
      //aux_in[i] = hits_[i];
    }

  auto cpu_input = cudautils::make_host_noncached_unique< HGCUncalibratedRecHit_GPU[] >(size_, cudaHostAllocWriteCombined);
  std::memcpy(cpu_input.get(), aux_in.data(), sizeof(HGCUncalibratedRecHit_GPU) * size_);
  auto d_oldhits = cudautils::make_device_unique<HGCUncalibratedRecHit_GPU[]>(size_, 0);
  auto d_newhits = cudautils::make_device_unique<HGCUncalibratedRecHit_GPU[]>(size_, 0);
  auto d_newhits_final = cudautils::make_device_unique<HGCRecHit_GPU[]>(size_, 0);
  auto h_newhits = cudautils::make_host_unique<HGCRecHit_GPU[]>(size_, 0);
  KernelManagerData<HGCUncalibratedRecHit_GPU, HGCRecHit_GPU> kmdata(cpu_input.get(), d_oldhits.get(), d_newhits.get(), d_newhits_final.get(), h_newhits.get());

  KernelManagerHGCalRecHit kernel_manager(kmdata);
  kernel_manager.run_kernels();

  //eventually this data will be produced in a 'convert_soa_data_to_collection()'
  //HGCRecHit_GPU* wrapper_out = kernel_manager.get_output();
  //edm::SortedCollection<T_OUT> coll = pointer_to_sorted_collection(wrapper_out, size_);
  //out_data_ = ...
}

template <class T_IN, class T_OUT>
edm::SortedCollection<T_OUT> HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::get_output_collection() 
{
  convert_soa_data_to_collection_<HGCRecHit>();
  //convert_collection_data_to_soa_<HGCRecHit>();
  return out_data_;
}

template <class T_IN, class T_OUT>
template <class U_IN>
void HeterogeneousProducerAcquireWrapper<T_IN, T_OUT>::convert_collection_data_to_soa_() 
{

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
  std::cout << "specialization!" << std::endl;
}

template class HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>;
template void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_soa_data_to_collection_<HGCRecHit>();
template void HeterogeneousProducerAcquireWrapper<HGCUncalibratedRecHit, HGCRecHit>::convert_collection_data_to_soa_<HGCUncalibratedRecHit>();
