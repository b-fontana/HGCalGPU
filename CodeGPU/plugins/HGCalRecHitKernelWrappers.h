#ifndef HGCalRecHitKernelWrappers_h
#define HGCalRecHitKernelWrappers_h

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "HeterogeneousCore/CUDACore/interface/CUDAScopedContext.h"
#include "HeterogeneousCore/CUDAUtilities/interface/device_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/host_noncached_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/host_unique_ptr.h"

#include "Utils.h"

template <typename T_IN, typename T_OUT>
  edm::SortedCollection<T_OUT> kernel_manager_wrapper(const edm::SortedCollection<T_IN>&, const CUDAScopedContextAcquire&);
							
#endif //HGCalRecHitKernelsWrappers_h
