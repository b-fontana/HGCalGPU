#include <cuda.h>
#include <cuda_runtime.h>

#include "KernelManager.h"
#include "HGCalRecHitKernelImpl.cuh"

dim3 nblocks_(1);
dim3 nthreads_(32);

KernelManagerHGCalRecHit::KernelManagerHGCalRecHit(const edm::SortedCollection<HGCUncalibratedRecHit>& oldhits):
  oldhits_collection_(oldhits)
{
  shits_ = oldhits_collection_.size();
  if(shits_ == 0)    
    throw cms::Exception("EmptyCollection") << "The collection is empty.";
  for(unsigned int i=0; i<shits_-1; ++i)
    assert(oldhits[i].id().det() == oldhits[i+1].id().det());
      
  DetId detid = oldhits[0].id();
  dtype_ = detid.det(); 
  if( (dtype_ != DetId::HGCalEE) and (dtype_ != DetId::HGCalHSi) and (dtype_ != DetId::HGCalHSc))
    throw cms::Exception("WrongDetectorType") << "The specified detector is wrong.";
    
  sbytes_ = shits_ * sizeof(oldhits_collection_[0]);

  cudaCheck( cudaMallocHost(&h_oldhits_, sbytes_) );
  cudaCheck( cudaMallocHost(&h_newhits_, sbytes_) );
  cudaCheck( cudaMalloc((void**)&d_oldhits_, sbytes_) );
  cudaCheck( cudaMalloc((void**)&d_newhits_, sbytes_) );
  cudaCheck( cudaMalloc((void**)&d_newhits_final_, sbytes_) );
  
  for(unsigned int j=0; j<shits_; ++j)
    h_oldhits_[j] = oldhits_collection_[j];
}

KernelManagerHGCalRecHit::~KernelManagerHGCalRecHit()
{
  cudaCheck( cudaFreeHost(h_oldhits_) );
  cudaCheck( cudaFreeHost(h_newhits_) );
  cudaCheck( cudaFree(d_oldhits_) );
  cudaCheck( cudaFree(d_newhits_) );
}

void KernelManagerHGCalRecHit::assign_and_transfer_to_device()
{
  cudaCheck( cudaMemcpyAsync(d_oldhits_, h_oldhits_, sbytes_, cudaMemcpyHostToDevice) );
  cudaCheck( cudaDeviceSynchronize() ); //needed because the copy is asynchronous
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::transfer_to_host_and_synchronize()
{
  cudaCheck( cudaMemcpyAsync(h_newhits_, d_newhits_final_, sbytes_, cudaMemcpyDeviceToHost) );
  cudaCheck( cudaDeviceSynchronize() );
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::reuse_device_pointers()
{
  std::swap(d_oldhits_, d_newhits_);
  cudaCheck( cudaDeviceSynchronize() );
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::run_kernels()
{
  assign_and_transfer_to_device();
  if(dtype_ == DetId::HGCalEE)
    {
      ee_step1_wrapper();
      reuse_device_pointers();
      to_rechit_wrapper();
    }
  else if(dtype_ == DetId::HGCalHSi)
    {
      hef_step1_wrapper();
      reuse_device_pointers();
      to_rechit_wrapper();
    }
  else
    {
      heb_step1_wrapper();
      reuse_device_pointers();  
      to_rechit_wrapper();
    }
  transfer_to_host_and_synchronize();
}

void KernelManagerHGCalRecHit::ee_step1_wrapper()
{
  ee_step1<<<nblocks_,nthreads_>>>(d_newhits_, d_oldhits_, shits_); 
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::hef_step1_wrapper()
{
  hef_step1<<<nblocks_,nthreads_>>>(d_newhits_, d_oldhits_, shits_); 
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::heb_step1_wrapper()
{  
  hef_step1<<<nblocks_,nthreads_>>>(d_newhits_, d_oldhits_, shits_); 
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::to_rechit_wrapper()
{
  to_rechit<<<nblocks_,nthreads_>>>(d_newhits_final_, d_oldhits_, shits_); 
  cudaCheck( cudaGetLastError() );
}

edm::SortedCollection<HGCRecHit> KernelManagerHGCalRecHit::get_new_collection()
{
  return pointer_to_sorted_collection<HGCRecHit>(h_newhits_, shits_);
}
