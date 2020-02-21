#include "Utils.h"

void print_to_histograms(HGCRecHitSoA *d, TH1F* &histo1, TH1F* &histo2, TH1F* &histo3, TH1F* &histo4, const unsigned int& nhits)
{
  for(unsigned int i = 0; i<nhits; ++i) 
    {
      histo1->Fill(d->energy[i]);
      histo2->Fill(d->time[i]);
      histo3->Fill(d->timeError[i]);
      histo4->Fill(d->son[i]);
    }
}

/*
template <typename T>
HGCRecHitCollection pointer_to_sorted_collection(T* ptr, const size_t& length)
{
  std::vector<T> v(ptr, ptr + length);
  HGCRecHitCollection coll(v);
  return coll;
}

template HGCRecHitCollection pointer_to_sorted_collection<HGCRecHit>(HGCRecHit*, const size_t&);
*/
