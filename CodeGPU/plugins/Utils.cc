#include "Utils.h"

/*
=======
>>>>>>> b5bfc7e2f47f926abb3dcd21cdf5e2094e53dd3f
template <typename T>
edm::SortedCollection<T> pointer_to_sorted_collection(T* ptr, const size_t& length)
{
  std::vector<T> v(ptr, ptr + length);
  edm::SortedCollection<T> coll(v);
  return coll;
}

<<<<<<< HEAD
template edm::SortedCollection<HGCRecHit_GPU> pointer_to_sorted_collection(HGCRecHit_GPU*, const size_t&);
*/
template <typename T>
HGCRecHitCollection pointer_to_sorted_collection(T* ptr, const size_t& length)
{
  std::vector<T> v(ptr, ptr + length);
  HGCRecHitCollection coll(v);
  return coll;
}

template HGCRecHitCollection pointer_to_sorted_collection<HGCRecHit>(HGCRecHit*, const size_t&);
