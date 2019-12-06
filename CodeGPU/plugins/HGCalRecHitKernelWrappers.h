#ifndef HGCalRecHitKernelWrappers_h
#define HGCalRecHitKernelWrappers_h

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "Utils.h"

template <typename T_IN, typename T_OUT>
  std::unique_ptr< edm::SortedCollection<T_OUT> > kernel_manager_wrapper(const edm::SortedCollection<T_IN>&,
									   const edm::EventSetup&);
							
#endif //HGCalRecHitKernelsWrappers_h
