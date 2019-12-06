#ifndef _UTILS_H_
#define _UTILS_H_

#include <vector>
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "HGCalRecHitKernelImpl.cuh"

enum detectortype { hgcee=0, hgchef=1, hgcheb=2, ntypes=3};

typedef edm::SortedCollection<HGCUncalibratedRecHit_GPU> HGCUncalibratedRecHitCollection_GPU;
typedef edm::SortedCollection<HGCRecHit_GPU> HGCRecHitCollection_GPU;

template <typename T>
HGCRecHitCollection pointer_to_sorted_collection(T*, const size_t&);

#endif //_UTILS_H_
