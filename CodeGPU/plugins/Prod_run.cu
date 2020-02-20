#include <cuda.h>
#include <cuda_runtime.h>
#include <inttypes.h>
#include "Prod_run.h"
#include "Kernel_impl.cuh"

dim3 nblocks_;
constexpr dim3 nthreads_(256); //some kernels will potentially not allocate shared memory properly with a lower number

Prod_run::Prod_run(const int& nhits, float *&h, float *&d1, float *&d2):
  nhits_(nhits), h_(h), d1_(d1), d2_(d2)
{
  nblocks_ = (nhits_ + nthreads_.x - 1) / nthreads_.x;
  nbytes = nhits_ * sizeof(float);
}

Prod_run::~Prod_run()
{
}

void Prod_run::assign_and_transfer_to_device_()
{
  cudaCheck( cudaMemcpyAsync(d1_, h_, 10, cudaMemcpyHostToDevice) );
  cudaCheck( cudaDeviceSynchronize() ); 
  cudaCheck( cudaGetLastError() );
}

void Prod_run::run_kernels()
{
  assign_and_transfer_to_device_();
  kernel<<<nblocks_, nthreads_>>>( d2_, d1_, nhits_ );
  after_kernel_();
}

void Prod_run::after_kernel_() {
  cudaCheck( cudaDeviceSynchronize() );
  cudaCheck( cudaGetLastError() );
}
