// -*- mode: c++; -*-

#ifndef XJETREC_UTIL_H_
#define XJETREC_UTIL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <cstdio>

namespace {

	void
	fprint_xmm_register(FILE *fp, const char *element_format,
						const int n)
	{
		char format[BUFSIZ];

		snprintf(format, BUFSIZ - 1,
				 "%%%%xmm%%d = [ %s | %s | %s | %s ]\n",
				 element_format, element_format, element_format,
				 element_format);

		unsigned int element_int[16] __attribute__ ((aligned(16)));
		float element_float[16] __attribute__ ((aligned(16)));

		switch(element_format[strlen(element_format) - 1]) {
		case 'd':
		case 'i':
		case 'o':
		case 'u':
		case 'x':
		case 'X':
			switch(n) {
			case 0:
				__asm__ __volatile__ ("movaps	%%xmm0, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm0");
				break;
			case 1:
				__asm__ __volatile__ ("movaps	%%xmm1, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm1");
				break;
			case 2:
				__asm__ __volatile__ ("movaps	%%xmm2, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm2");
				break;
			case 3:
				__asm__ __volatile__ ("movaps	%%xmm3, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm3");
				break;
			case 4:
				__asm__ __volatile__ ("movaps	%%xmm4, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm4");
				break;
			case 5:
				__asm__ __volatile__ ("movaps	%%xmm5, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm5");
				break;
			case 6:
				__asm__ __volatile__ ("movaps	%%xmm6, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm6");
				break;
			case 7:
				__asm__ __volatile__ ("movaps	%%xmm7, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm7");
				break;
#ifdef __x86_64__
			case 8:
				__asm__ __volatile__ ("movaps	%%xmm8, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm8");
				break;
			case 9:
				__asm__ __volatile__ ("movaps	%%xmm9, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm9");
				break;
			case 10:
				__asm__ __volatile__ ("movaps	%%xmm10, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm10");
				break;
			case 11:
				__asm__ __volatile__ ("movaps	%%xmm11, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm11");
				break;
			case 12:
				__asm__ __volatile__ ("movaps	%%xmm12, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm12");
				break;
			case 13:
				__asm__ __volatile__ ("movaps	%%xmm13, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm13");
				break;
			case 14:
				__asm__ __volatile__ ("movaps	%%xmm14, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm14");
				break;
			case 15:
				__asm__ __volatile__ ("movaps	%%xmm15, %0"
					: "=m" (*element_int) : "m" (*element_int)
					: "%xmm15");
				break;
#endif // __x86_64__
			}
			fprintf(fp, format, n,
					element_int[3], element_int[2],
					element_int[1], element_int[0]);
			break;
		case 'a':
		case 'A':
		case 'e':
		case 'E':
		case 'f':
		case 'F':
		case 'g':
		case 'G':
			switch(n) {
			case 0:
				__asm__ __volatile__ ("movaps	%%xmm0, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm0");
				break;
			case 1:
				__asm__ __volatile__ ("movaps	%%xmm1, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm1");
				break;
			case 2:
				__asm__ __volatile__ ("movaps	%%xmm2, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm2");
				break;
			case 3:
				__asm__ __volatile__ ("movaps	%%xmm3, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm3");
				break;
			case 4:
				__asm__ __volatile__ ("movaps	%%xmm4, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm4");
				break;
			case 5:
				__asm__ __volatile__ ("movaps	%%xmm5, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm5");
				break;
			case 6:
				__asm__ __volatile__ ("movaps	%%xmm6, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm6");
				break;
			case 7:
				__asm__ __volatile__ ("movaps	%%xmm7, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm7");
				break;
#ifdef __x86_64__
			case 8:
				__asm__ __volatile__ ("movaps	%%xmm8, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm8");
				break;
			case 9:
				__asm__ __volatile__ ("movaps	%%xmm9, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm9");
				break;
			case 10:
				__asm__ __volatile__ ("movaps	%%xmm10, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm10");
				break;
			case 11:
				__asm__ __volatile__ ("movaps	%%xmm11, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm11");
				break;
			case 12:
				__asm__ __volatile__ ("movaps	%%xmm12, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm12");
				break;
			case 13:
				__asm__ __volatile__ ("movaps	%%xmm13, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm13");
				break;
			case 14:
				__asm__ __volatile__ ("movaps	%%xmm14, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm14");
				break;
			case 15:
				__asm__ __volatile__ ("movaps	%%xmm15, %0"
					: "=m" (*element_float) : "m" (*element_float)
					: "%xmm15");
				break;
#endif // __x86_64__
			}
			fprintf(fp, format, n,
					element_float[3], element_float[2],
					element_float[1], element_float[0]);
			break;
		}
	}
}

#endif // XJETREC_UTIL_H_
