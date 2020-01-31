#include <algorithm>
#include <cuda.h>
#include <cuda_runtime.h>

#include "KernelManager.h"
#include "HGCalRecHitKernelImpl.cuh"

dim3 nblocks_;
constexpr dim3 nthreads_(256);

KernelManagerHGCalRecHit::KernelManagerHGCalRecHit(KernelModifiableData<HGCUncalibratedRecHitSoA, HGCRecHitSoA> data, const DetId::Detector& dtype):
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

void KernelManagerHGCalRecHit::assign_and_transfer_to_device(const KernelConstantData<HGCeeUncalibratedRecHitConstantData>& h_kcdata, const KernelConstantData<HGCeeUncalibratedRecHitConstantData>& d_kcdata)
{
  cudaCheck( cudaMemcpyAsync( d_kcdata.data.hgcEE_fCPerMIP_, h_kcdata.data.hgcEE_fCPerMIP_, h_kcdata.data.nbytes, cudaMemcpyHostToDevice) );
  cudaCheck( cudaDeviceSynchronize() ); //needed because the copy is asynchronous
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::assign_and_transfer_to_device(const KernelConstantData<HGChefUncalibratedRecHitConstantData>& h_kcdata, const KernelConstantData<HGChefUncalibratedRecHitConstantData>& d_kcdata)
{
  cudaCheck( cudaMemcpyAsync( d_kcdata.data.hgcHEF_fCPerMIP_, h_kcdata.data.hgcHEF_fCPerMIP_, h_kcdata.vdata.nbytes, cudaMemcpyHostToDevice) );
  cudaCheck( cudaDeviceSynchronize() ); //needed because the copy is asynchronous
  cudaCheck( cudaGetLastError() );
}

void KernelManagerHGCalRecHit::assign_and_transfer_to_device(const KernelConstantData<HGChebUncalibratedRecHitConstantData>& h_kcdata, const KernelConstantData<HGChebUncalibratedRecHitConstantData>& d_kcdata)
{
  cudaCheck( cudaMemcpyAsync( d_kcdata.data.hgcHEB_fCPerMIP_, h_kcdata.data.hgcHEB_fCPerMIP_, h_kcdata.data.nbytes, cudaMemcpyHostToDevice) );
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

void KernelManagerHGCalRecHit::run_kernels(const KernelConstantData<HGCeeUncalibratedRecHitConstantData>& h_kcdata, const KernelConstantData<HGCeeUncalibratedRecHitConstantData>& d_kcdata)
{
  printf("%d blocks being launched with %d threads (%d in total).\n", nblocks_.x, nthreads_.x, nblocks_.x*nthreads_.x);
  assign_and_transfer_to_device(h_kcdata, d_kcdata);
  assign_and_transfer_to_device();

  printf("Running ee kernel with: %zu hits.\n", data_.nhits);
  ee_step1<<<nblocks_,nthreads_>>>( *(data_.d_2), *(data_.d_1), d_kcdata.data, data_.nhits);
  after_kernel();

  reuse_device_pointers();

  ee_to_rechit<<<nblocks_,nthreads_>>>( *(data_.d_out), *(data_.d_1), d_kcdata.data, data_.nhits );
  after_kernel();

  transfer_to_host_and_synchronize();
}

void KernelManagerHGCalRecHit::run_kernels(const KernelConstantData<HGChefUncalibratedRecHitConstantData>& h_kcdata, const KernelConstantData<HGChefUncalibratedRecHitConstantData>& d_kcdata)
{
  printf("%d blocks being launched with %d threads (%d in total).\n", nblocks_.x, nthreads_.x, nblocks_.x*nthreads_.x);
  assign_and_transfer_to_device(h_kcdata, d_kcdata);
  assign_and_transfer_to_device();

  printf("Running ee kernel with: %zu hits.\n", data_.nhits);
  hef_step1<<<nblocks_,nthreads_>>>( *(data_.d_2), *(data_.d_1), d_kcdata.data, data_.nhits);
  after_kernel();

  reuse_device_pointers();

  hef_to_rechit<<<nblocks_,nthreads_>>>( *(data_.d_out), *(data_.d_1), d_kcdata.data, data_.nhits );
  after_kernel();

  transfer_to_host_and_synchronize();
}

void KernelManagerHGCalRecHit::run_kernels(const KernelConstantData<HGChebUncalibratedRecHitConstantData>& h_kcdata, const KernelConstantData<HGChebUncalibratedRecHitConstantData>& d_kcdata)
{
  printf("%d blocks being launched with %d threads (%d in total).\n", nblocks_.x, nthreads_.x, nblocks_.x*nthreads_.x);
  assign_and_transfer_to_device(h_kcdata, d_kcdata);
  assign_and_transfer_to_device();

  printf("Running ee kernel with: %zu hits.\n", data_.nhits);
  heb_step1<<<nblocks_,nthreads_>>>( *(data_.d_2), *(data_.d_1), d_kcdata.data, data_.nhits);
  after_kernel();

  reuse_device_pointers();

  heb_to_rechit<<<nblocks_,nthreads_>>>( *(data_.d_out), *(data_.d_1), d_kcdata.data, data_.nhits );
  after_kernel();

  transfer_to_host_and_synchronize();
}

void KernelManagerHGCalRecHit::after_kernel() {
  cudaCheck( cudaDeviceSynchronize() );
  cudaCheck( cudaGetLastError() );
}

HGCRecHitSoA* KernelManagerHGCalRecHit::get_output()
{
  return data_.h_out;
}
