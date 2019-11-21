#include <cstdio>
#include <iostream>
#include <cuda_runtime.h>

#include "HGCalRecHitKernel.cuh"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

__global__
void scrambler_kernel(const char * message, size_t length)
{
  //printf("blockIdx.x, threadIdx.x: %d, %d\n", blockIdx.x, threadIdx.x);
  for (size_t i = blockDim.x * blockIdx.x + threadIdx.x; i < length; i += blockDim.x * gridDim.x)
    printf("%c", message[i]);
}

__host__
void scrambler_wrapper(const char * message, size_t length) {
  scrambler_kernel<<<16,1>>>(message, length);
  cudaCheck(cudaDeviceSynchronize());
  cudaCheck(cudaGetLastError());
  std::cout << std::endl;

  scrambler_kernel<<<4,4>>>(message, length);
  cudaCheck(cudaDeviceSynchronize());
  cudaCheck(cudaGetLastError());
  std::cout << std::endl;

  scrambler_kernel<<<1,16>>>(message, length);
  cudaCheck(cudaDeviceSynchronize());
  cudaCheck(cudaGetLastError());
  std::cout << std::endl;
}

__global__
void HGCeeRecHitKernel(HGCUncalibratedRecHit* oldhits, 
		       HGCUncalibratedRecHit* newhits, size_t length)
{
  for (size_t i = blockDim.x * blockIdx.x + threadIdx.x; i < length; i += blockDim.x * gridDim.x)
    {
      newhits[i] = oldhits[i];
    }
}

__host__
std::unique_ptr<HGCeeUncalibratedRecHitCollection> HGCeeRecHitKernel_wrapper(HGCeeUncalibratedRecHitCollection eeHits, 
									      const std::string& collection_name)
{
  auto eeRecHits = std::make_unique<HGCeeUncalibratedRecHitCollection>();

  HGCUncalibratedRecHit h_oldhits[eeHits.size()], h_newhits[eeHits.size()];

  size_t sz = eeHits.size() * sizeof(eeHits[0]);
  for(size_t i=0; i<eeHits.size(); ++i)
    h_oldhits[i] = eeHits[i];

  HGCUncalibratedRecHit *d_oldhits, *d_newhits;
  cudaCheck( cudaMalloc((void**)&d_oldhits, sz) );
  cudaCheck( cudaMemcpy(d_oldhits, h_oldhits, sz, cudaMemcpyHostToDevice) );
  cudaCheck( cudaMalloc((void**)&d_newhits, sz) );

  std::cout << sz << ", " << eeHits.size() << std::endl;
  HGCeeRecHitKernel<<<1,32>>>(d_oldhits, d_newhits, eeHits.size());

  cudaCheck(cudaDeviceSynchronize());
  cudaCheck(cudaGetLastError());
  std::cout << sz << ", " << eeHits.size() << std::endl;
  cudaCheck( cudaMemcpy(h_newhits, d_newhits, sz, cudaMemcpyDeviceToHost) );
  std::cout << sz << ", " << eeHits.size() << std::endl;
  std::cout << "DATA:" << std::endl;
  std::cout << h_newhits[0].amplitude() << std::endl;
  std::vector<HGCUncalibratedRecHit> newhits_vec(h_newhits, h_newhits + eeHits.size());
  HGCeeUncalibratedRecHitCollection newhits_coll(newhits_vec);
  auto newhits = std::make_unique<HGCeeUncalibratedRecHitCollection>(newhits_coll);

  cudaCheck( cudaFree(d_oldhits) );
  cudaCheck( cudaFree(d_newhits) );
  return newhits;
}

