// -*- mode: c++; -*-

#ifndef XJETREC_MEM_H_
#define XJETREC_MEM_H_

#include <cstddef>
#include <cstdlib>
#include <cerrno>

namespace {

	/**
	 * Aligned memory allocation
	 *
	 * With the following two exceptions, this implementation conforms
	 * to IEEE Std 1003.1d-1999:
	 *
	 * The non-standard function builtin_aligned_free() deallocates
	 * memory that has previously been allocated by builtin_memalign()
	 * (IEEE Std 1003.1d-1999, line 170-171).
	 *
	 * Alignment may be not a power of two multiple of sizeof(void *),
	 * and EINVAL is never returned (IEEE Std 1003.1d-1999, line
	 * 182-183).
	 *
	 * @param[out] memptr pointer to the allocated memory
	 * @param[in] alignment a non-zero alignment boundary that may be
	 * not a power of two multiple of sizeof(void *)
	 * @param[in] size allocated bytes
	 * @see IEEE Std 1003.1d-1999, section 19.2.2 (p. 67f)
	 */
	int builtin_memalign(
		void **memptr, size_t alignment, size_t size)
	{
		*memptr = malloc(sizeof(size_t) + alignment + size);

		// Required by IEEE Std 1003.1d-1999, line 182-183
		if(*memptr == NULL) {
			// IEEE Std 1003.1d-1999 technically permits an undefined
			// *memptr upon returning ENOMEM.
			return ENOMEM;
		}

		const size_t remainder =
			(*reinterpret_cast<size_t *>(memptr) + sizeof(size_t)) %
			alignment;
		size_t offset = sizeof(size_t) +
			(remainder == 0 ? 0 : alignment - remainder);

		*reinterpret_cast<size_t *>(memptr) += offset;
		*(reinterpret_cast<size_t *>(*memptr) - 1) = offset;

		// Required by IEEE Std 1003.1d-1999, line 176-177
		return 0;
	}

	/**
	 * Deallocates the aligned space pointed to ptr
	 *
	 * With the following exception, this implementation conforms
	 * to ISO/IEC 9899:1999:
	 *
	 * If the argument does not match a pointer earlier returned by
	 * builtin_memalign, or if the space has been deallocated by a
	 * call to builtin_aligned_free, the behavior is undefined.
	 *
	 * @param[in] ptr aligned space pointer
	 * @see ISO/IEC 9899:1999, section 7.20.3.2
	 */
	void builtin_aligned_free(void *ptr)
	{
		// Required by ISO/IEC 9899:1999, section 7.20.3.2
		if(ptr == NULL) {
			return;
		}

		const size_t offset =
			*(reinterpret_cast<size_t *>(ptr) - 1);

		free(reinterpret_cast<unsigned char *>(ptr) - offset);
	}

}

#endif // XJETREC_MEM_H_
