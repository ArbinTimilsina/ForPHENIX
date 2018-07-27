// This file is part of PPMd project
// Written and distributed to public domain by Dmitry Shkarin 1997,
// 1999-2001, 2006
// Contents: main routine
// Comments: system & compiler dependent file

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <ppmd/PPMd.h>

namespace ppmd {

#define USAGE_STR \
	"Usage: PPMd <e|d> [switches] <FileName[s] | Wildcard[s]>\n"
	static const char *pFName;
	static uint32_t StartFilePosition;
	static bool EncodeFlag;
	static clock_t StartClock;
	static struct ARC_INFO
	{
		// FileLength & CRC? Hmm, maybe in another times...
		uint32_t signature, attrib;
		uint16_t info, FNLen, time, date;
	} ai;

#include <sys/stat.h>
#include <time.h>
#include <utime.h>
#include <dirent.h>
#include <unistd.h>
#include <fnmatch.h>

	inline void EnvSetNormAttr(const char *FName)
	{
		chmod(FName, S_IWUSR|S_IRUSR);
	}
	inline int EnvGetCh(void)
	{
		return getchar();
	}
	inline void EnvGetCWD(char *CurDir)
	{
		getcwd(CurDir, 260);
	}
	inline void EnvSetDateTimeAttr(const char *WrkStr)
	{
		struct utimbuf t;
		t.actime = t.modtime = (ai.date << 16) + ai.time;
		utime(WrkStr, &t);
		chmod(WrkStr, ai.attrib);
	}
	struct ENV_FIND_RESULT
	{
		dirent de;
		struct stat st;
		const char *getFName() const
		{
			return de.d_name;
		}
		void copyDateTimeAttr() const
		{
			ai.attrib = st.st_mode;
			ai.time = st.st_mtime & 0xFFFF;
			ai.date = st.st_mtime >> 16;
		}
	};
	struct ENV_FILE_FINDER
	{
		const char *pPattern;
		DIR *dir;
		dirent *de;
		struct stat st;
		ENV_FIND_RESULT getResult()
		{
			ENV_FIND_RESULT Rslt;
			Rslt.de = *de;
			Rslt.st = st;
			return Rslt;
		}
		bool isFileValid()
		{
			return (fnmatch(pPattern, de->d_name,
							FNM_NOESCAPE) == 0 &&
					stat(de->d_name, &st) == 0 &&
					(st.st_mode & S_IRUSR) != 0 && st.st_nlink == 1);
		}
		bool findFirst(const char *Pattern)
		{
			pPattern = Pattern;
			return ((dir = opendir(".")) &&
					(de = readdir(dir)) != NULL);
		}
		bool findNext()
		{
			return ((de = readdir(dir)) != NULL);
		}
		void findStop()
		{
			closedir(dir);
		}
	};

	static const char *const MTxt[] = {
		"Can`t open file %s\n",
		"read/write error for files %s/%s\n", "Out of memory!\n",
		"User break\n", "unknown command: %s\n",
		"unknown switch: %s\n",
		"designed and written by Dmitry Shkarin "
		"<dmitry.shkarin@mtu-net.ru>\n"
		USAGE_STR
		"Switches (for encoding only):\n"
		"\t-d	  - delete file[s] after processing, default: "
		"disabled\n"
		"\t-fName - set output file name to Name\n"
		"\t-mN	  - use N MB memory - [1,256], default: %d\n"
		"\t-oN	  - set model order to N - [2,%d], default: %d\n"
		"\t-rN	  - set method of model restoration at memory "
		"insufficiency:\n"
		"\t\t-r0 - restart model from scratch (default)\n"
		"\t\t-r1 - cut off model (slow)\n"
	};

	void PrintInfo(_PPMD_FILE *DecodedFile, _PPMD_FILE *EncodedFile)
	{
		char WrkStr[320];
		uint32_t NDec = ftell(DecodedFile);
		NDec += (NDec == 0);
		uint32_t NEnc = ftell(EncodedFile) - StartFilePosition;
		uint32_t n1 = (8U * NEnc) / NDec;
		uint32_t n2 =
			(100U * (8U * NEnc - NDec * n1) + NDec / 2U) / NDec;
		if (n2 == 100) {
			n1++;
			n2 = 0;
		}
		int RunTime = ((clock() - StartClock) << 10) /
			int(CLOCKS_PER_SEC);
		uint32_t Speed = NDec / (RunTime + (RunTime == 0));
		uint32_t UsedMemory = GetUsedMemory() >> 10;
		uint32_t m1 = UsedMemory >> 10;
		uint32_t m2 =
			(10U * (UsedMemory - (m1 << 10)) + (1 << 9)) >> 10;
		if (m2 == 10) {
			m1++;
			m2 = 0;
		}
		if(!EncodeFlag)
			SWAP(NDec, NEnc);
		sprintf(WrkStr, "%14s:%7d >%7d, %1d.%02d bpb, used:%3d.%1d "
				"MB, speed: %d KB/sec",
				pFName, (int)NDec, (int)NEnc, (int)n1, (int)n2,
				(int)m1, (int)m2, (int)Speed);
		printf("%-79.79s\r", WrkStr);
		fflush(stdout);
	}

	static char *ChangeExtRare(const char *In, char *Out,
							   const char *Ext)
	{
		char *RetVal = Out;
		const char *p = strrchr(In, '.');
		if(!p || strrchr(In, '/') > p)
			p = In + strlen(In);
		do {
			*Out++ = *In++;
		}
		while(In != p);
		*Out++ = '.';
		while((*Out++ = *Ext++) != 0)
			;
		return RetVal;
	}
	inline bool RemoveFile(const char *FName)
	{
		EnvSetNormAttr(FName);
		return (remove(FName) == 0);
	}
	static bool TestAccessRare(const char *FName)
	{
		static bool YesToAll = false;
		FILE *fp = fopen(FName, "rb");
		if(!fp)
			return true;
		fclose(fp);
		if(YesToAll)
			return RemoveFile(FName);
		printf("%s already exists, overwrite?: <Y>es, <N>o, <A>ll, "
			   "<Q>uit?", FName);
		for(;;)
			switch(toupper(EnvGetCh())) {
			case 'A':
				YesToAll = true;
			case '\r':
			case 'Y':
				return RemoveFile(FName);
			case 0x1B:
			case 'Q':
				printf(MTxt[3]);
				exit(-1);
			case 'N':
				return false;
			}
	}
	static FILE *FOpen(const char *FName, const char *mode)
	{
		FILE *fp = fopen(FName, mode);
		if(!fp) {
			printf(MTxt[0], FName);
			exit(-1);
		}
		setvbuf(fp, NULL, _IOFBF, 64 * 1024);
		return fp;
	}
	inline void PrepareCoding(int SASize, FILE *fp)
	{
		if(!StartSubAllocator(SASize)) {
			printf(MTxt[2]);
			exit(-1);
		}
		StartClock = clock();
		StartFilePosition = ftell(fp);
	}
	inline void EncodeFile(const ENV_FIND_RESULT& efr, int MaxOrder,
						   int SASize, bool CutOff,
						   const char *ArcName)
	{
		char WrkStr[260];
		strcpy(WrkStr, ArcName);
		if(!WrkStr[0] &&
		   !TestAccessRare(ChangeExtRare(efr.getFName(), WrkStr,
										 "pmd")))
			return;
		FILE *fpIn = FOpen(efr.getFName(), "rb");
		FILE *fpOut = FOpen(WrkStr, "a+b");
		pFName = strrchr(efr.getFName(), '/');
		pFName = (pFName) ? (pFName + 1) : (efr.getFName());
		efr.copyDateTimeAttr();
		ai.signature = PPMdSignature;
		ai.FNLen = strlen(pFName) + (CutOff << 14);
		ai.info = (MaxOrder - 1) | ((SASize - 1) << 4) |
			((PROG_VAR - 'A') << 12);
		fwrite(&ai, sizeof(ai), 1, fpOut);
		fwrite(pFName, ai.FNLen & 0x1FF, 1, fpOut);
		PrepareCoding(SASize, fpOut);
		EncodeFile(fpOut, fpIn, MaxOrder, CutOff);
		putchar('\n');
		if(ferror(fpOut) || ferror(fpIn)) {
			printf(MTxt[1], efr.getFName(), WrkStr);
			exit(-1);
		}
		fclose(fpIn);
		fclose(fpOut);
	}
	inline bool DecodeOneFile(FILE *fpIn)
	{
		char WrkStr[260];
		int MaxOrder, SASize;
		bool CutOff;
		if(!fread(&ai, sizeof(ai), 1, fpIn))
			return false;
		CutOff = ai.FNLen >> 14;
		ai.FNLen = CLAMP(int(ai.FNLen & 0x1FF), 1, 260 - 1);
		fread(WrkStr, ai.FNLen, 1, fpIn);
		WrkStr[ai.FNLen] = 0;
		if(!TestAccessRare(WrkStr))
			return false;
		FILE *fpOut = FOpen(pFName = WrkStr, "wb");
		MaxOrder = (ai.info & 0x0F) + 1;
		SASize = ((ai.info >> 4) & 0xFF) + 1;
		uint32_t Variant = (ai.info >> 12) + 'A';
		if(ai.signature != PPMdSignature || Variant != PROG_VAR) {
			printf(MTxt[0], WrkStr);
			exit(-1);
		}
		PrepareCoding(SASize, fpIn);
		DecodeFile(fpOut, fpIn, MaxOrder, CutOff);
		putchar('\n');
		if(ferror(fpOut) || ferror(fpIn) || feof(fpIn)) {
			printf(MTxt[1], WrkStr, WrkStr);
			exit(-1);
		}
		fclose(fpOut);
		EnvSetDateTimeAttr(WrkStr);
		return true;
	}
	inline void DecodeFile(const ENV_FIND_RESULT& efr)
	{
		FILE *fpIn = FOpen(efr.getFName(), "rb");
		while(DecodeOneFile(fpIn))
			;
		fclose(fpIn);
	}
	inline void TestArchive(char *ArcName, const char *Pattern)
	{
		if(!Pattern[0]) {
			char CurDir[260];
			EnvGetCWD(CurDir);
			const char *p = strrchr(CurDir, '/');
			p = (p && strlen(p + 1)) ? (p + 1) : ("PPMdFile");
			ChangeExtRare(p, ArcName, "pmd");
		}
		else
			strcpy(ArcName, Pattern);
		FILE *fp = fopen(ArcName, "rb");
		if(fp) {
			if(!fread(&ai, sizeof(ai), 1, fp) ||
			   ai.signature != PPMdSignature ||
			   (ai.info >> 12) + 'A' != PROG_VAR) {
				printf(MTxt[0], ArcName);
				exit(-1);
			}
			fclose(fp);
		}
	}
	struct FILE_LIST_NODE
	{
		FILE_LIST_NODE *next;
		ENV_FIND_RESULT efr;
		FILE_LIST_NODE(const ENV_FIND_RESULT &Data,
					   FILE_LIST_NODE **PrevNode)
		{
			efr = Data;
			next = *PrevNode;
			*PrevNode = this;
		}
		void destroy(FILE_LIST_NODE **PrevNode)
		{
			*PrevNode = next;
			delete this;
		}
	};
}

int main(int argc, char *argv[])
{
	assert(TestCompilation());
	char ArcName[260];
	bool DeleteFile = false, CutOff = false;
	int i, MaxOrder = 4, SASize = 10;
	printf("Fast PPMII compressor for textual data, variant %c, "
		   __DATE__ "\n", ppmd::PROG_VAR);
	if (argc < 3) {
		printf(ppmd::MTxt[6], SASize, ppmd::MAX_O, MaxOrder);
		return -1;
	}
	switch(toupper(argv[1][0])) {
	case 'E':
		ppmd::EncodeFlag = true;
		break;
	case 'D':
		ppmd::EncodeFlag = false;
		break;
	default:
		printf(ppmd::MTxt[4], argv[1]);
		return -1;
	}
	for(ArcName[0] = 0, i = 2;
		i < argc && (argv[i][0] == '-' || argv[i][0] == '/'); i++)
		switch(toupper(argv[i][1])) {
		case 'D':
			DeleteFile = true;
			break;
		case 'F':
			ppmd::TestArchive(ArcName, argv[i] + 2);
			break;
		case 'M':
			SASize = ppmd::CLAMP(atoi(argv[i] + 2), 1, 256);
			break;
		case 'O':
			MaxOrder = ppmd::CLAMP(atoi(argv[i] + 2), 2,
								   int(ppmd::MAX_O));
			break;
		case 'R':
			CutOff = ppmd::CLAMP(atoi(argv[i] + 2), 0, 1);
			break;
		default:
			printf(ppmd::MTxt[5], argv[i]);
			return -1;
		}
	ppmd::FILE_LIST_NODE *pNode;
	ppmd::FILE_LIST_NODE *pFirstNode = NULL;
	ppmd::FILE_LIST_NODE **ppNode = &pFirstNode;

	for(ppmd::ENV_FILE_FINDER eff; i < argc; i++) {
		if(eff.findFirst(argv[i]))
			do {
				if(eff.isFileValid()) {
					pNode = new ppmd::FILE_LIST_NODE(eff.getResult(), ppNode);
					if(!pNode) {
						printf(ppmd::MTxt[2]);
						return -1;
					}
					ppNode = &(pNode->next);
				}
			}
			while(eff.findNext());
		eff.findStop();
	}
	while((pNode = pFirstNode) != NULL) {
		ppmd::ENV_FIND_RESULT &efr = pNode->efr;
		if(ppmd::EncodeFlag)
			ppmd::EncodeFile(efr, MaxOrder, SASize, CutOff, ArcName);
		else
			ppmd::DecodeFile(efr);
		if(DeleteFile)
			remove(efr.getFName());
		pNode->destroy(&pFirstNode);
	}
	ppmd::StopSubAllocator();
	return 0;
}
