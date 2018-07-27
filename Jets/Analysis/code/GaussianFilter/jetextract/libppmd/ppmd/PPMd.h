// This file is part of PPMd project
// Written and distributed to public domain by Dmitry Shkarin 1997,
// 1999-2001, 2006
// Contents: interface to encoding/decoding routines
// Comments: this file can be used as an interface to PPMd module
// (consisting of Model.cpp) from external program

#ifndef _PPMD_H_
#define _PPMD_H_

#include <ppmd/PPMdType.h>

#ifdef __cplusplus
extern "C"
{
    namespace ppmd
    {
#endif // __cplusplus
		bool StartSubAllocator(uint32_t SubAllocatorSize);
		void StopSubAllocator(void);	// It can be called once
		uint32_t GetUsedMemory(void);	// for information only

		// (MaxOrder == 1) parameter value has special meaning, it
		// does not restart model and can be used for solid mode
		// archives; Call sequence:
		//     StartSubAllocator(SubAllocatorSize);
		//     EncodeFile(SolidArcFile, File1, MaxOrder, true);
		//     EncodeFile(SolidArcFile, File2,        1, true);
		//     ...
		//     EncodeFile(SolidArcFile, FileN,        1, true);
		//     StopSubAllocator(void);
		void EncodeFile(_PPMD_FILE *EncodedFile,
						_PPMD_FILE *DecodedFile,
						int MaxOrder, bool CutOff);
		void DecodeFile(_PPMD_FILE *DecodedFile,
						_PPMD_FILE *EncodedFile,
						int MaxOrder, bool CutOff);

		// imported function
		void PrintInfo(_PPMD_FILE *DecodedFile,
					   _PPMD_FILE *EncodedFile);
#ifdef __cplusplus
    }
}
#endif // __cplusplus
#endif // _PPMD_H_
