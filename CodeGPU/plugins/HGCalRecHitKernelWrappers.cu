#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>
#include <type_traits>
#include <cuda_runtime.h>

#include "HGCalRecHitKernelWrappers.h"

template <typename T>
__host__
std::unique_ptr<T> hgc_rechit_kernel_wrapper(T hits, Detector det)
{
  if (!std::is_same<T, HGCUncalibratedRecHitCollection>::value)
    {
      std::cout << "The hgc_rechit_kernel_wrapper template does not currently support this type." << std::endl;
      std::exit(1);
    }

  size_t shits = hits.size();
  if(shits==0)
    return std::make_unique<T>(); //owns nothing
  size_t sz = shits * sizeof(hits[0]);

  HGCUncalibratedRecHit *h_oldhits, *h_newhits;
  cudaCheck( cudaMallocHost(&h_oldhits, sz) );
  cudaCheck( cudaMallocHost(&h_newhits, sz) );

  for(size_t i=0; i<shits; ++i)
    h_oldhits[i] = hits[i];

  HGCUncalibratedRecHit *d_oldhits, *d_newhits;
  cudaCheck( cudaMalloc((void**)&d_oldhits, sz) );
  cudaCheck( cudaMemcpy(d_oldhits, h_oldhits, sz, cudaMemcpyHostToDevice) );
  cudaCheck( cudaMalloc((void**)&d_newhits, sz) );

  dim3 dimGrid(1);
  dim3 dimBlock(32);
  if (det == Detector::EE)
    HGCRecHitKernel<<<dimGrid,dimBlock>>>(d_oldhits, d_newhits, shits);
  else if (det == Detector::HEF)
    HGCRecHitKernel<<<dimGrid,dimBlock>>>(d_oldhits, d_newhits, shits);
  else if (det == Detector::HEB)
    HGCRecHitKernel<<<dimGrid,dimBlock>>>(d_oldhits, d_newhits, shits);
  else {
    std::cout << "Error: HGCalRecHitKernelWrappers.cu" << std::endl;
    std::exit(1);
  }

  cudaCheck( cudaDeviceSynchronize() );
  cudaCheck( cudaGetLastError() );
  cudaCheck( cudaMemcpy(h_newhits, d_newhits, sz, cudaMemcpyDeviceToHost) );

  std::unique_ptr<T> newhits = std::make_unique<T>( pointer_to_sorted_collection<HGCUncalibratedRecHit>(h_newhits, shits) );

  cudaCheck( cudaFreeHost(h_oldhits) );
  cudaCheck( cudaFreeHost(h_newhits) );
  cudaCheck( cudaFree(d_oldhits) );
  cudaCheck( cudaFree(d_newhits) );
 
  return newhits;
}

template std::unique_ptr<HGCUncalibratedRecHitCollection> hgc_rechit_kernel_wrapper<HGCUncalibratedRecHitCollection>(HGCUncalibratedRecHitCollection, Detector);
