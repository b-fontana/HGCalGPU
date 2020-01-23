#include <algorithm>
#include <cuda.h>
#include <cuda_runtime.h>

#include "KernelManager.h"
#include "HGCalRecHitKernelImpl.cuh"

dim3 nblocks_;
constexpr dim3 nthreads_(256);

KernelManagerHGCalRecHit::KernelManagerHGCalRecHit(KernelManagerData<HGCUncalibratedRecHitSoA, HGCRecHitSoA> data, const DetId::Detector& dtype):
  data_(data), dtype_(dtype)
{
  nblocks_ = (data_.nhits + nthreads_.x - 1) / nthreads_.x; 
}

KernelManagerHGCalRecHit::~KernelManagerHGCalRecHit()
{
}

void KernelManagerHGCalRecHit::assign_and_transfer_to_device()
{
  cudaCheck( cudaMemcpyAsync((data_.d_1)->amplitude, (data_.h_in)->amplitude, (data_.d_1)->nbytes, cudaMemcpyHostToDevice) );
  cudaCheck( cudaDeviceSynchronize() ); //needed because the copy is asynchronous
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::transfer_to_host_and_synchronize()
{
  cudaCheck( cudaMemcpyAsync((data_.h_out)->energy, (data_.d_out)->energy, (data_.d_out)->nbytes, cudaMemcpyDeviceToHost) );
  cudaCheck( cudaDeviceSynchronize() );
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::reuse_device_pointers()
{
  std::swap(data_.d_1, data_.d_2); 
  cudaCheck( cudaDeviceSynchronize() );
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::run_kernels()
{
  assign_and_transfer_to_device();

  if(dtype_ == DetId::HGCalEE)
    {
      std::cout << "to_rechit() start" << std::endl;
      //ee_step1_wrapper();
      //reuse_device_pointers();
      //to_rechit_wrapper();
      std::cout << "to_rechit() end" << std::endl;
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
  ee_step1<<<nblocks_,nthreads_>>>(data_.d_2, data_.d_1, data_.nhits); 
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::hef_step1_wrapper()
{
  hef_step1<<<nblocks_,nthreads_>>>(data_.d_2, data_.d_1, data_.nhits); 
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::heb_step1_wrapper()
{  
  hef_step1<<<nblocks_,nthreads_>>>(data_.d_2, data_.d_1, data_.nhits); 
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::to_rechit_wrapper()
{
  to_rechit<<<nblocks_,nthreads_>>>(data_.d_out, data_.d_1, data_.nhits); 
  cudaCheck( cudaGetLastError() );
}

HGCRecHitSoA* KernelManagerHGCalRecHit::get_output()
{
  return data_.h_out;
}
