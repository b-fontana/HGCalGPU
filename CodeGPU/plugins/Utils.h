#ifndef _UTILS_H_
#define _UTILS_H_

#include <vector>
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "Types.h"

enum detectortype { hgcee=0, hgchef=1, hgcheb=2, ntypes=3};

typedef edm::SortedCollection<HGCUncalibratedRecHitSoA> HGCUncalibratedRecHitCollectionSoA;
typedef edm::SortedCollection<HGCRecHitSoA> HGCRecHitCollectionSoA;

/*
template <typename T>
edm::SortedCollection<T> pointer_to_sorted_collection(T*, const size_t&);
*/

#endif //_UTILS_H_
