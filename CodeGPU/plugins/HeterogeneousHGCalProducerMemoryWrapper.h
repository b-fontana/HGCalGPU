#ifndef _HETEROGENEOUSHGCALPRODUCERMEMORYWRAPPER_H_
#define _HETEROGENEOUSHGCALPRODUCERMEMORYWRAPPER_H_

#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>
#include <type_traits>
#include <numeric>
#include <cuda_runtime.h>

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "HeterogeneousCore/CUDACore/interface/CUDAScopedContext.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAContextState.h"
#include "HeterogeneousCore/CUDAServices/interface/CUDAService.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HeterogeneousCore/CUDAUtilities/interface/device_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/host_noncached_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/host_unique_ptr.h"

#include "KernelManager.h"
#include "Utils.h"
#include "Types.h"


namespace memory {
  namespace allocation {
    /*
    namespace {
      std::tuple<LENGTHSIZE, LENGTHSIZE, LENGTHSIZE, LENGTHSIZE> get_memory_sizes_(const std::vector<LENGTHSIZE>&, const LENGTHSIZE&, const LENGTHSIZE&, const LENGTHSIZE&);
    }
    */
    void host(KernelConstantData<HGCeeUncalibratedRecHitConstantData>*, cudautils::host::noncached::unique_ptr<double[]>&);
    void host(KernelConstantData<HGChefUncalibratedRecHitConstantData>*, cudautils::host::noncached::unique_ptr<double[]>&);
    void host(KernelConstantData<HGChebUncalibratedRecHitConstantData>*, cudautils::host::noncached::unique_ptr<double[]>&);
    void host(const int&, HGCUncalibratedRecHitSoA*, cudautils::host::noncached::unique_ptr<float[]>&);
    void host(const int&, HGCRecHitSoA*&, cudautils::host::unique_ptr<float[]>&);
    void device(KernelConstantData<HGCeeUncalibratedRecHitConstantData>*, cudautils::device::unique_ptr<double[]>&);
    void device(KernelConstantData<HGChefUncalibratedRecHitConstantData>*, cudautils::device::unique_ptr<double[]>&);
    void device(KernelConstantData<HGChebUncalibratedRecHitConstantData>*, cudautils::device::unique_ptr<double[]>&);
    void device(const int&, HGCUncalibratedRecHitSoA*, HGCUncalibratedRecHitSoA*, HGCRecHitSoA*, cudautils::device::unique_ptr<float[]>&);
  }
}
							
#endif // _HETEROGENEOUSHGCALPRODUCERMEMORYWRAPPER_H_
