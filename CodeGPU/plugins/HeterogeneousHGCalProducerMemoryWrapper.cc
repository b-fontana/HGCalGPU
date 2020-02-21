#include "HeterogeneousHGCalProducerMemoryWrapper.h"
#include "Types.h"

namespace memory {
  namespace allocation {    
    namespace {
      //returns total number of bytes, number of 'double' elements and number of 'float' elements
      std::tuple<LENGTHSIZE, LENGTHSIZE, LENGTHSIZE, LENGTHSIZE> get_memory_sizes_(const std::vector<LENGTHSIZE>& fixed_sizes, const LENGTHSIZE& ndoubles, const LENGTHSIZE& nfloats, const LENGTHSIZE& nints)
      {
	const LENGTHSIZE size1 = sizeof(double);
	const LENGTHSIZE size2 = sizeof(float);
	const LENGTHSIZE size3 = sizeof(int);
	LENGTHSIZE nelements1_tot = std::accumulate( fixed_sizes.begin(), fixed_sizes.begin() + ndoubles, 0);
	LENGTHSIZE nelements2_tot = std::accumulate( fixed_sizes.begin() + ndoubles, fixed_sizes.begin() + ndoubles + nfloats, 0);
	LENGTHSIZE nelements3_tot = std::accumulate( fixed_sizes.begin() + ndoubles + nfloats, fixed_sizes.end(), 0);
	assert( fixed_sizes.begin() + ndoubles + nfloats + nints == fixed_sizes.end() );
	LENGTHSIZE size_tot = nelements1_tot*size1+nelements2_tot*size2+nelements3_tot*size3;
	return std::make_tuple(size_tot, nelements1_tot, nelements2_tot, nelements3_tot);
      }
    }

    void device(KernelConstantData<HGCeeUncalibratedRecHitConstantData> *kcdata, cudautils::device::unique_ptr<double[]>& mem) {
      const std::vector<LENGTHSIZE> nelements = {kcdata->data.s_hgcEE_fCPerMIP_, kcdata->data.s_hgcEE_cce_, kcdata->data.s_hgcEE_noise_fC_, kcdata->data.s_rcorr_, kcdata->data.s_weights_, kcdata->data.s_waferTypeL_};
      auto memsizes = get_memory_sizes_(nelements, 5, 0, 1);

      mem = cudautils::make_device_unique<double[]>(std::get<0>(memsizes), 0);

      kcdata->data.hgcEE_fCPerMIP_ = mem.get();
      kcdata->data.hgcEE_cce_      = kcdata->data.hgcEE_fCPerMIP_ + nelements[0];
      kcdata->data.hgcEE_noise_fC_ = kcdata->data.hgcEE_cce_ + nelements[1];
      kcdata->data.rcorr_          = kcdata->data.hgcEE_noise_fC_ + nelements[2];
      kcdata->data.weights_        = kcdata->data.rcorr_ + nelements[3];
      kcdata->data.waferTypeL_     = reinterpret_cast<int*>(kcdata->data.weights_ + nelements[4]);
      kcdata->data.nbytes = std::get<0>(memsizes);
      kcdata->data.ndelem = std::get<1>(memsizes) + 2;
      kcdata->data.nfelem = std::get<2>(memsizes) + 4;
      kcdata->data.nielem = std::get<3>(memsizes) + 0;
      kcdata->data.nuelem = 2;
      kcdata->data.nbelem = 1;
    }

    void device(KernelConstantData<HGChefUncalibratedRecHitConstantData> *kcdata, cudautils::device::unique_ptr<double[]>& mem) {
      const std::vector<LENGTHSIZE> nelements = {kcdata->data.s_hgcHEF_fCPerMIP_, kcdata->data.s_hgcHEF_cce_, kcdata->data.s_hgcHEF_noise_fC_, kcdata->data.s_rcorr_, kcdata->data.s_weights_, kcdata->data.s_waferTypeL_};
      auto memsizes = get_memory_sizes_(nelements, 5, 0, 1);

      mem = cudautils::make_device_unique<double[]>(std::get<0>(memsizes), 0);

      kcdata->data.hgcHEF_fCPerMIP_ = mem.get();
      kcdata->data.hgcHEF_cce_      = kcdata->data.hgcHEF_fCPerMIP_ + nelements[0];
      kcdata->data.hgcHEF_noise_fC_ = kcdata->data.hgcHEF_cce_ + nelements[1];
      kcdata->data.rcorr_           = kcdata->data.hgcHEF_noise_fC_ + nelements[2];
      kcdata->data.weights_         = kcdata->data.rcorr_ + nelements[3];
      kcdata->data.waferTypeL_      = reinterpret_cast<int*>(kcdata->data.weights_ + nelements[4]);
      kcdata->data.nbytes = std::get<0>(memsizes);
      kcdata->data.ndelem = std::get<1>(memsizes) + 2;
      kcdata->data.nfelem = std::get<2>(memsizes) + 4;
      kcdata->data.nielem = std::get<3>(memsizes) + 0;
      kcdata->data.nuelem = 3;
      kcdata->data.nbelem = 1;
    }

    void device(KernelConstantData<HGChebUncalibratedRecHitConstantData> *kcdata, cudautils::device::unique_ptr<double[]>& mem) {
      const std::vector<LENGTHSIZE> nelements = {kcdata->data.s_weights_};
      auto memsizes = get_memory_sizes_(nelements, 1, 0, 0);

      mem = cudautils::make_device_unique<double[]>(std::get<0>(memsizes), 0);

      kcdata->data.weights_         = mem.get();
      kcdata->data.nbytes = std::get<0>(memsizes);
      kcdata->data.ndelem = std::get<1>(memsizes) + 3;
      kcdata->data.nfelem = std::get<2>(memsizes) + 0;
      kcdata->data.nielem = std::get<3>(memsizes) + 0;
      kcdata->data.nuelem = 3;
      kcdata->data.nbelem = 1;
    }

    void device(const int& nhits, HGCUncalibratedRecHitSoA* soa1, HGCUncalibratedRecHitSoA* soa2, HGCRecHitSoA* soa3, cudautils::device::unique_ptr<float[]>& mem)
    {
      std::vector<LENGTHSIZE> sizes = {6*sizeof(float), 3*sizeof(uint32_t),                     //soa1
				       6*sizeof(float), 3*sizeof(uint32_t),                     //soa2
				       3*sizeof(float), 2*sizeof(uint32_t), 1*sizeof(uint8_t)}; //soa3
      LENGTHSIZE size_tot = std::accumulate( sizes.begin(), sizes.end(), 0);

      mem = cudautils::make_device_unique<float[]>(nhits * size_tot, 0);

      soa1->amplitude     = mem.get();
      soa1->pedestal      = soa1->amplitude    + nhits;
      soa1->jitter        = soa1->pedestal     + nhits;
      soa1->chi2          = soa1->jitter       + nhits;
      soa1->OOTamplitude  = soa1->chi2         + nhits;
      soa1->OOTchi2       = soa1->OOTamplitude + nhits;
      soa1->flags         = reinterpret_cast<uint32_t*>(soa1->OOTchi2 + nhits);
      soa1->aux           = soa1->flags        + nhits;
      soa1->id            = soa1->aux          + nhits;

      soa2->amplitude     = reinterpret_cast<float*>(soa1->id + nhits);
      soa2->pedestal      = soa2->amplitude    + nhits;
      soa2->jitter        = soa2->pedestal     + nhits;
      soa2->chi2          = soa2->jitter       + nhits;
      soa2->OOTamplitude  = soa2->chi2         + nhits;
      soa2->OOTchi2       = soa2->OOTamplitude + nhits;
      soa2->flags         = reinterpret_cast<uint32_t*>(soa2->OOTchi2 + nhits);
      soa2->aux           = soa2->flags        + nhits;
      soa2->id            = soa2->aux          + nhits;
  
      soa3->energy        = reinterpret_cast<float*>(soa2->id + nhits);
      soa3->time          = soa3->energy       + nhits;
      soa3->timeError     = soa3->time         + nhits;
      soa3->id            = reinterpret_cast<uint32_t*>(soa2->id + nhits);
      soa3->flagBits      = soa3->id           + nhits;
      soa3->son           = reinterpret_cast<uint8_t*>(soa3->flagBits + nhits);

      soa1->nbytes = std::accumulate(sizes.begin(), sizes.begin()+2, 0);
      soa2->nbytes = std::accumulate(sizes.begin()+2, sizes.begin()+4, 0);
      soa3->nbytes = std::accumulate(sizes.begin()+4, sizes.end(), 0);
    }

    void host(KernelConstantData<HGCeeUncalibratedRecHitConstantData>* kcdata, cudautils::host::noncached::unique_ptr<double[]>& mem)
    {
      const std::vector<LENGTHSIZE> nelements = {kcdata->data.s_hgcEE_fCPerMIP_, kcdata->data.s_hgcEE_cce_, kcdata->data.s_hgcEE_noise_fC_, kcdata->data.s_rcorr_, kcdata->data.s_weights_, kcdata->data.s_waferTypeL_};
      auto memsizes = get_memory_sizes_(nelements, 5, 0, 1);

      mem = cudautils::make_host_noncached_unique<double[]>(std::get<0>(memsizes), 0);

      kcdata->data.hgcEE_fCPerMIP_ = mem.get();
      kcdata->data.hgcEE_cce_      = kcdata->data.hgcEE_fCPerMIP_ + nelements[0];
      kcdata->data.hgcEE_noise_fC_ = kcdata->data.hgcEE_cce_ + nelements[1];
      kcdata->data.rcorr_          = kcdata->data.hgcEE_noise_fC_ + nelements[2];
      kcdata->data.weights_        = kcdata->data.rcorr_ + nelements[3];
      kcdata->data.waferTypeL_     = reinterpret_cast<int*>(kcdata->data.weights_ + nelements[4]);
      kcdata->data.nbytes = std::get<0>(memsizes);
      kcdata->data.ndelem = std::get<1>(memsizes) + 2;
      kcdata->data.nfelem = std::get<2>(memsizes) + 0;
      kcdata->data.nielem = std::get<3>(memsizes) + 0;
      kcdata->data.nuelem = 2;
      kcdata->data.nbelem = 1;
    }

    void host(KernelConstantData<HGChefUncalibratedRecHitConstantData>* kcdata, cudautils::host::noncached::unique_ptr<double[]>& mem)
    {
      const std::vector<LENGTHSIZE> nelements = {kcdata->data.s_hgcHEF_fCPerMIP_, kcdata->data.s_hgcHEF_cce_, kcdata->data.s_hgcHEF_noise_fC_, kcdata->data.s_rcorr_, kcdata->data.s_weights_, kcdata->data.s_waferTypeL_};
      auto memsizes = get_memory_sizes_(nelements, 5, 0, 1);

      mem = cudautils::make_host_noncached_unique<double[]>(std::get<0>(memsizes), 0);

      kcdata->data.hgcHEF_fCPerMIP_ = mem.get();
      kcdata->data.hgcHEF_cce_      = kcdata->data.hgcHEF_fCPerMIP_ + nelements[0];
      kcdata->data.hgcHEF_noise_fC_ = kcdata->data.hgcHEF_cce_ + nelements[1];
      kcdata->data.rcorr_           = kcdata->data.hgcHEF_noise_fC_ + nelements[2];
      kcdata->data.weights_         = kcdata->data.rcorr_ + nelements[3];
      kcdata->data.waferTypeL_      = reinterpret_cast<int*>(kcdata->data.weights_ + nelements[4]);
      kcdata->data.nbytes = std::get<0>(memsizes);
      kcdata->data.ndelem = std::get<1>(memsizes) + 2;
      kcdata->data.nfelem = std::get<2>(memsizes) + 0;
      kcdata->data.nielem = std::get<3>(memsizes) + 0;
      kcdata->data.nuelem = 3;
      kcdata->data.nbelem = 1;
    }

    void host(KernelConstantData<HGChebUncalibratedRecHitConstantData>* kcdata, cudautils::host::noncached::unique_ptr<double[]>& mem)
    {
      const std::vector<LENGTHSIZE> nelements = {kcdata->data.s_weights_};
      auto memsizes = get_memory_sizes_(nelements, 1, 0, 0);

      mem = cudautils::make_host_noncached_unique<double[]>(std::get<0>(memsizes), 0);

      kcdata->data.weights_         = mem.get();
      kcdata->data.nbytes = std::get<0>(memsizes);
      kcdata->data.ndelem = std::get<1>(memsizes) + 3;
      kcdata->data.nfelem = std::get<2>(memsizes) + 0;
      kcdata->data.nielem = std::get<3>(memsizes) + 0;
      kcdata->data.nuelem = 3;
      kcdata->data.nbelem = 1;
    }

    void host(const int& nhits, HGCUncalibratedRecHitSoA* soa, cudautils::host::noncached::unique_ptr<float[]>& mem)
    {
      LENGTHSIZE size1 = (LENGTHSIZE)6*sizeof(float);
      LENGTHSIZE size2 = (LENGTHSIZE)3*sizeof(uint32_t);
      LENGTHSIZE size_tot = size1 + size2;
      mem = cudautils::make_host_noncached_unique<float[]>(nhits * size_tot, 0);

      soa->amplitude     = (float*)(mem.get());
      soa->pedestal      = (float*)(soa->amplitude    + nhits);
      soa->jitter        = (float*)(soa->pedestal     + nhits);
      soa->chi2          = (float*)(soa->jitter       + nhits);
      soa->OOTamplitude  = (float*)(soa->chi2         + nhits);
      soa->OOTchi2       = (float*)(soa->OOTamplitude + nhits);
      soa->flags         = (uint32_t*)(soa->OOTchi2   + nhits);
      soa->aux           = (uint32_t*)(soa->flags     + nhits);
      soa->id            = (uint32_t*)(soa->aux       + nhits);
      soa->nbytes = size_tot;
    }

    void host(const int& nhits, HGCRecHitSoA*& soa, cudautils::host::unique_ptr<float[]>& mem)
    {
      LENGTHSIZE size1 = (LENGTHSIZE)(3*sizeof(float));
      LENGTHSIZE size2 = (LENGTHSIZE)(2*sizeof(uint32_t));
      LENGTHSIZE size3 = (LENGTHSIZE)(1*sizeof(uint8_t));
      LENGTHSIZE size_tot = size1 + size2 + size3;
      mem = cudautils::make_host_unique<float[]>(nhits*size_tot, 0);

      soa->energy     = (float*)(mem.get());
      soa->time       = (float*)(soa->energy     + nhits);
      soa->timeError  = (float*)(soa->time       + nhits);
      soa->id         = (uint32_t*)(soa->time    + nhits);
      soa->flagBits   = (uint32_t*)(soa->id      + nhits);
      soa->son        = (uint8_t*)(soa->flagBits + nhits);
      soa->nbytes = size_tot;
    }
  }
}
