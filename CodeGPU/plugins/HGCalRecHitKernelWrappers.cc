#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>
#include <type_traits>
#include <cuda_runtime.h>

#include "HGCalRecHitKernelWrappers.h"
#include "KernelManager.h"
#include "HGCalRecHitKernelImpl.cuh"

/*
I am currently mixing GPU and non-GPU types 
Eventually all types will be given by T_IN and T_OUT
*/
template <typename T_IN, typename T_OUT>
edm::SortedCollection<T_OUT> kernel_manager_wrapper(const edm::SortedCollection<T_IN>& hits, const CUDAScopedContextAcquire& ctx)
{
  if (!std::is_same<T_IN, HGCUncalibratedRecHit>::value)
    throw cms::Exception("WrongTemplateType") << "The hgc_rechit_kernel_wrapper template does not support this type.";
  
  size_t size = hits.size();
  if (size == 0)
    return edm::SortedCollection<T_OUT>();

  for(unsigned int i=0; i<size-1; ++i)
    assert(hits[i].id().det() == hits[i+1].id().det());

  DetId::Detector dtype_ = hits[0].id().det(); 
  if( (dtype_ != DetId::HGCalEE) and (dtype_ != DetId::HGCalHSi) and (dtype_ != DetId::HGCalHSc))
    throw cms::Exception("WrongDetectorType") << "The specified detector is wrong.";

  std::vector<HGCUncalibratedRecHit_GPU> aux_in(size);
  for(uint i=0; i<size; ++i)
    {
      HGCUncalibratedRecHit_GPU hit;
      aux_in[i] = hit;
      //aux_in[i] = hits[i];
    }

  auto cpu_input = cudautils::make_host_noncached_unique< HGCUncalibratedRecHit_GPU[] >(size, cudaHostAllocWriteCombined);
  std::memcpy(cpu_input.get(), aux_in.data(), sizeof(HGCUncalibratedRecHit_GPU) * size);
  auto d_oldhits = cudautils::make_device_unique<HGCUncalibratedRecHit_GPU[]>(size, ctx.stream());
  auto d_newhits = cudautils::make_device_unique<HGCUncalibratedRecHit_GPU[]>(size, ctx.stream());
  auto d_newhits_final = cudautils::make_device_unique<HGCRecHit_GPU[]>(size, ctx.stream());
  auto h_newhits = cudautils::make_host_unique<HGCRecHit_GPU[]>(size, ctx.stream());

  KernelManagerData<HGCUncalibratedRecHit_GPU, HGCRecHit_GPU> data(cpu_input.get(), d_oldhits.get(), d_newhits.get(), d_newhits_final.get(), h_newhits.get());
  KernelManagerHGCalRecHit kernel_manager(data);
  kernel_manager.run_kernels();

  //HGCRecHit_GPU* wrapper_out = kernel_manager.get_output();
  //edm::SortedCollection<T_OUT> coll = pointer_to_sorted_collection(wrapper_out, size);
  edm::SortedCollection<HGCRecHit> coll;

  return coll;
}

template edm::SortedCollection<HGCRecHit> kernel_manager_wrapper<HGCUncalibratedRecHit, HGCRecHit>( const edm::SortedCollection<HGCUncalibratedRecHit>&, const CUDAScopedContextAcquire&);
