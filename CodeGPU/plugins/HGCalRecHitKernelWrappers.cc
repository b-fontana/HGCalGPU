#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>
#include <type_traits>
#include <cuda_runtime.h>

#include "HGCalRecHitKernelWrappers.h"
#include "KernelManager.h"

template <typename T_IN, typename T_OUT>
std::unique_ptr< edm::SortedCollection<T_OUT> > kernel_manager_wrapper(const edm::SortedCollection<T_IN>& hits, const edm::EventSetup& es)
{
  if (!std::is_same<T_IN, HGCUncalibratedRecHit>::value)
    throw cms::Exception("WrongTemplateType") << "The hgc_rechit_kernel_wrapper template does not support this type.";

  if (hits.size() == 0)
    return std::make_unique< edm::SortedCollection<T_OUT> >();

  KernelManagerHGCalRecHit kernel_manager(hits);
  kernel_manager.run_kernels();

  edm::SortedCollection<T_OUT> coll = kernel_manager.get_new_collection();
  auto newhits = std::make_unique< edm::SortedCollection<T_OUT> >(coll);
  return newhits;
}

template std::unique_ptr< edm::SortedCollection<HGCRecHit> > kernel_manager_wrapper<HGCUncalibratedRecHit, HGCRecHit>( const edm::SortedCollection<HGCUncalibratedRecHit>&, const edm::EventSetup&);
