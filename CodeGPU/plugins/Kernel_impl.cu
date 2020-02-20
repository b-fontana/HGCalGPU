#include <cuda.h>
#include <cuda_runtime.h>

__global__
void kernel(float *d2, float *d1, int length)
{
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  for (int i = tid; i < length; i += blockDim.x * gridDim.x)
    d2[i] = d1[i] * 2;
}
