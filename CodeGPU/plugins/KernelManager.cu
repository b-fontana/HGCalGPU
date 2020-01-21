#include <algorithm>
#include <cuda.h>
#include <cuda_runtime.h>

#include "KernelManager.h"
#include "HGCalRecHitKernelImpl.cuh"

dim3 nblocks_(1);
dim3 nthreads_(32);

KernelManagerHGCalRecHit::KernelManagerHGCalRecHit(KernelManagerData<HGCUncalibratedRecHitSoA, HGCRecHitSoA> data, const DetId::Detector& dtype):
  data_(data), dtype_(dtype)
{
  sbytes_ = shits_ * sizeof(data_.h_in_[0]);
}

KernelManagerHGCalRecHit::~KernelManagerHGCalRecHit()
{
}

void KernelManagerHGCalRecHit::assign_and_transfer_to_device()
{
  cudaCheck( cudaMemcpyAsync(data_.d_1_, data_.h_in_, sbytes_, cudaMemcpyHostToDevice) );
  cudaCheck( cudaDeviceSynchronize() ); //needed because the copy is asynchronous
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::transfer_to_host_and_synchronize()
{
  cudaCheck( cudaMemcpyAsync(data_.h_out_, data_.d_out_, sbytes_, cudaMemcpyDeviceToHost) );
  cudaCheck( cudaDeviceSynchronize() );
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::reuse_device_pointers()
{
  //std::swap(data_.d_1_, data_.d_2_); 
  cudaCheck( cudaDeviceSynchronize() );
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::run_kernels()
{
  assign_and_transfer_to_device();
 
  /*
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
  */

  transfer_to_host_and_synchronize();
}

void KernelManagerHGCalRecHit::ee_step1_wrapper()
{
  ee_step1<<<nblocks_,nthreads_>>>(data_.d_2_, data_.d_1_, shits_); 
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::hef_step1_wrapper()
{
  hef_step1<<<nblocks_,nthreads_>>>(data_.d_2_, data_.d_1_, shits_); 
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::heb_step1_wrapper()
{  
  hef_step1<<<nblocks_,nthreads_>>>(data_.d_2_, data_.d_1_, shits_); 
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::to_rechit_wrapper()
{
  to_rechit<<<nblocks_,nthreads_>>>(data_.d_out_, data_.d_1_, shits_); 
  cudaCheck( cudaGetLastError() );
}

HGCRecHitSoA* KernelManagerHGCalRecHit::get_output()
{
  return data_.h_out_;
}
