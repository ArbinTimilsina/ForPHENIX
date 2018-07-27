// -*- mode: c; -*-

#ifndef XJETREC_CUDA_H_
#define XJETREC_CUDA_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>
#include <vector_types.h>

extern void cufft_multiply_scale_32_256(float2 a[], const float2 b[],
										int size, float scale);
extern void cuda_device_init(void);
extern void cuda_malloc(void **device_ptr, size_t size);
extern void cuda_memcpy(void *destination, const void *source,
						size_t count, enum cudaMemcpyKind kind);
extern void cuda_free(void *device_ptr);
extern void cufft_plan_1d(cufftHandle *plan, int nx, cufftType type,
						  int batch);
extern void cufft_plan_2d(cufftHandle *plan, int nx, int ny,
						  cufftType type);
extern void cufft_plan_3d(cufftHandle *plan, int nx, int ny, int nz,
						  cufftType type);
extern void cufft_destroy(cufftHandle plan);
extern void cufft_exec_c2c(cufftHandle plan,
						   cufftComplex *input_data,
						   cufftComplex *output_data, int direction);
extern void cufft_exec_r2c(cufftHandle plan, cufftReal *input_data,
						   cufftComplex *output_data);
extern void cufft_exec_c2r(cufftHandle plan,
						   cufftComplex *input_data,
						   cufftReal *output_data);

#endif // XJETREC_CUDA_H_
