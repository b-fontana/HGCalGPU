#ifndef _PRODRUN_H
#define _PRODRUN_H

#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCompat.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "Kernel_impl.cuh"
#include "Types.h"

#include <vector>
#include <algorithm> //std::swap  
#include <variant>
#include <cuda.h>
#include <cuda_runtime.h>

#ifdef __CUDA_ARCH__
extern __constant__ uint32_t calo_rechit_masks[];
#endif

class Prod_run {
 public:
  Prod_run(const int&, float*&, float*&, float*&);
  ~Prod_run();
  void run_kernels();

 private:
  void after_kernel_();
  void assign_and_transfer_to_device_();

  int nbytes;
  int nhits_;
  float *h_;
  float *d1_;
  float *d2_;
};

#endif //_PRODRUN_H_
