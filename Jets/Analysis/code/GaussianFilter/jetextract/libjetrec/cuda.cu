// -*- mode: c; -*-

#include <stdio.h>
#include <jetrec/cuda.h>

#define CUDA_CHECK_SUCESS(result) \
	if(result != CUDA_SUCCESS) { \
        fprintf(stderr, "%s:%d: CUDA driver error %x\n", \
				__FILE__, __LINE__, result); \
		switch(result) { \
		case CUDA_ERROR_INVALID_VALUE: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_INVALID_VALUE"); \
			break; \
		case CUDA_ERROR_OUT_OF_MEMORY: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_OUT_OF_MEMORY"); \
			break; \
		case CUDA_ERROR_NOT_INITIALIZED: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_NOT_INITIALIZED"); \
			break; \
		case CUDA_ERROR_NO_DEVICE: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_NO_DEVICE"); \
			break; \
		case CUDA_ERROR_INVALID_DEVICE: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_INVALID_DEVICE"); \
			break; \
		case CUDA_ERROR_INVALID_IMAGE: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_INVALID_IMAGE"); \
			break; \
		case CUDA_ERROR_INVALID_CONTEXT: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_INVALID_CONTEXT"); \
			break; \
		case CUDA_ERROR_CONTEXT_ALREADY_CURRENT: \
			fprintf(stderr, " (%s)", \
					"CUDA_ERROR_CONTEXT_ALREADY_CURRENT"); \
			break; \
		case CUDA_ERROR_MAP_FAILED: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_MAP_FAILED"); \
			break; \
		case CUDA_ERROR_UNMAP_FAILED: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_UNMAP_FAILED"); \
			break; \
		case CUDA_ERROR_ARRAY_IS_MAPPED: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_ARRAY_IS_MAPPED"); \
			break; \
		case CUDA_ERROR_ALREADY_MAPPED: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_ALREADY_MAPPED"); \
			break; \
		case CUDA_ERROR_NO_BINARY_FOR_GPU: \
			fprintf(stderr, " (%s)", \
					"CUDA_ERROR_NO_BINARY_FOR_GPU"); \
			break; \
		case CUDA_ERROR_ALREADY_ACQUIRED: \
			fprintf(stderr, " (%s)", \
					"CUDA_ERROR_ALREADY_ACQUIRED"); \
			break; \
		case CUDA_ERROR_NOT_MAPPED: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_NOT_MAPPED"); \
			break; \
		case CUDA_ERROR_INVALID_SOURCE: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_INVALID_SOURCE"); \
			break; \
		case CUDA_ERROR_FILE_NOT_FOUND: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_FILE_NOT_FOUND"); \
			break; \
		case CUDA_ERROR_INVALID_HANDLE: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_INVALID_HANDLE"); \
			break; \
		case CUDA_ERROR_NOT_FOUND: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_NOT_FOUND"); \
			break; \
		case CUDA_ERROR_LAUNCH_FAILED: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_LAUNCH_FAILED"); \
			break; \
		case CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES: \
			fprintf(stderr, " (%s)", \
					"CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES"); \
			break; \
		case CUDA_ERROR_LAUNCH_TIMEOUT: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_LAUNCH_TIMEOUT"); \
			break; \
		case CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING: \
			fprintf(stderr, " (%s)", \
					"CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING"); \
			break; \
		case CUDA_ERROR_UNKNOWN: \
			fprintf(stderr, " (%s)", "CUDA_ERROR_UNKNOWN"); \
			break; \
		}; \
		fprintf(stderr, "\n"); \
        exit(EXIT_FAILURE); \
	}

#define CUFFT_CHECK_SUCESS(result) \
	if(result != CUFFT_SUCCESS) { \
        fprintf(stderr, "%s:%d: CUFFT driver error %x", \
				__FILE__, __LINE__, result); \
		switch(result) { \
		case CUFFT_INVALID_PLAN: \
			fprintf(stderr, " (%s)", "CUFFT_INVALID_PLAN"); \
			break; \
		case CUFFT_ALLOC_FAILED: \
			fprintf(stderr, " (%s)", "CUFFT_ALLOC_FAILED"); \
			break; \
		case CUFFT_INVALID_TYPE: \
			fprintf(stderr, " (%s)", "CUFFT_INVALID_TYPE"); \
			break; \
		case CUFFT_INVALID_VALUE: \
			fprintf(stderr, " (%s)", "CUFFT_INVALID_VALUE"); \
			break; \
		case CUFFT_INTERNAL_ERROR: \
			fprintf(stderr, " (%s)", "CUFFT_INTERNAL_ERROR"); \
			break; \
		case CUFFT_EXEC_FAILED: \
			fprintf(stderr, " (%s)", "CUFFT_EXEC_FAILED"); \
			break; \
		case CUFFT_SETUP_FAILED: \
			fprintf(stderr, " (%s)", "CUFFT_SETUP_FAILED"); \
			break; \
		case CUFFT_INVALID_SIZE: \
			fprintf(stderr, " (%s)", "CUFFT_INVALID_SIZE"); \
			break; \
		}; \
		fprintf(stderr, "\n"); \
        exit(EXIT_FAILURE); \
	}

static __device__ __host__ inline float2
complex_scale(const float2 a, const float s)
{
    float2 c = {a.x * s, a.y * s};

    return c;
}

static __device__ __host__ inline float2
complex_multiply(const float2 a, const float2 b)
{
    float2 c = {a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x};

    return c;
}

static __global__ void
multiply_scale(float2 a[], const float2 b[], int size, float scale)
{
    const int num_threads = blockDim.x * gridDim.x;
    const int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    for(int i = thread_id; i < size; i += num_threads)
        a[i] = complex_scale(complex_multiply(a[i], b[i]), scale);     
}

void cufft_multiply_scale_32_256(float2 a[], const float2 b[],
								 int size, float scale)
{
	multiply_scale<<<32, 256>>>(a, b, size, scale);

	const cudaError_t result = cudaThreadSynchronize();

	CUDA_CHECK_SUCESS(result);
}

void cuda_device_init(void)
{
	cudaError_t result;
    int device_count;
    int dev;

	result = cudaGetDeviceCount(&device_count);
	CUDA_CHECK_SUCESS(result);
    if(device_count == 0) {
        fprintf(stderr, "%s:%d: no device found\n",
				__FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    for(dev = 0; dev < device_count; dev++) {
        cudaDeviceProp device_prop;

        result = cudaGetDeviceProperties(&device_prop, dev);
		CUDA_CHECK_SUCESS(result);
        if(device_prop.major >= 1)
            break;
    }
    if(dev == device_count) {
        fprintf(stderr, "%s:%d: no device supporting CUDA found\n",
				__FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    else {
        result = cudaSetDevice(dev);
		CUDA_CHECK_SUCESS(result);
	}
}

void cuda_malloc(void **device_ptr, size_t size)
{
	const cudaError_t result = cudaMalloc(device_ptr, size);

	CUDA_CHECK_SUCESS(result);
}

void cuda_memcpy(void *destination, const void *source,
				 size_t count, enum cudaMemcpyKind kind)
{
	const cudaError_t result =
		cudaMemcpy(destination, source, count, kind);

	CUDA_CHECK_SUCESS(result);
}

void cuda_free(void *device_ptr)
{
	const cudaError_t result = cudaFree(device_ptr);

	CUDA_CHECK_SUCESS(result);
}

void cufft_plan_1d(cufftHandle *plan, int nx, cufftType type,
				   int batch)
{
	const cufftResult result = cufftPlan1d(plan, nx, type, batch);

	CUFFT_CHECK_SUCESS(result);
}

void cufft_plan_2d(cufftHandle *plan, int nx, int ny,
				   cufftType type)
{
	const cufftResult result = cufftPlan2d(plan, nx, ny, type);

	CUFFT_CHECK_SUCESS(result);
}

void cufft_plan_3d(cufftHandle *plan, int nx, int ny, int nz,
				   cufftType type)
{
	const cufftResult result = cufftPlan3d(plan, nx, ny, nz, type);

	CUFFT_CHECK_SUCESS(result);
}

void cufft_destroy(cufftHandle plan)
{
	const cufftResult result = cufftDestroy(plan);

	CUFFT_CHECK_SUCESS(result);
}

void cufft_exec_c2c(cufftHandle plan, cufftComplex *input_data,
					cufftComplex *output_data, int direction)
{
	const cufftResult result =
		cufftExecC2C(plan, input_data, output_data, direction);

	CUFFT_CHECK_SUCESS(result);
}

void cufft_exec_r2c(cufftHandle plan, cufftReal *input_data,
					cufftComplex *output_data)
{
	const cufftResult result =
		cufftExecR2C(plan, input_data, output_data);

	CUFFT_CHECK_SUCESS(result);
}

void cufft_exec_c2r(cufftHandle plan, cufftComplex *input_data,
					cufftReal *output_data)
{
	const cufftResult result =
		cufftExecC2R(plan, input_data, output_data);

	CUFFT_CHECK_SUCESS(result);
}

