// -*- mode: c++; -*-

// This file is part of PPMd project
// Written and distributed to public domain by Dmitry Shkarin 1997,
// 1999-2001, 2006
// Contents: memory allocation routines

namespace ppmd {

	enum {
		UNIT_SIZE = 12, N1 = 4, N2 = 4, N3 = 4,
		N4 = (128 + 3 - 1 * N1 - 2 * N2 - 3 * N3) / 4,
		N_INDICES = N1 + N2 + N3 + N4
	};

	inline void PrefetchData(void *Addr)
	{
#ifdef _USE_PREFETCHING
		// FIXME: Assembly NetBurst prefetching needed
		// uint8_t PrefetchByte = *(volatile uint8_t *)Addr;
#endif // _USE_PREFETCHING
	}
	// constants
	static uint8_t Indx2Units[N_INDICES], Units2Indx[128];
	static __thread uint32_t GlueCount, GlueCount1,
		SubAllocatorSize = 0;
	static __thread uint8_t *HeapStart, *pText, *UnitsStart;
	static __thread uint8_t *LoUnit, *HiUnit, *AuxUnit;

#if defined(_32_NORMAL) || defined(_64_EXOTIC)
	inline uint32_t Ptr2Indx(void *p)
	{
		return (uint32_t)p;
	}
	inline void *Indx2Ptr(uint32_t indx)
	{
		return (void *)indx;
	}
#else // defined(_32_NORMAL) || defined(_64_EXOTIC)
	static __thread uint8_t *HeapNull;
	inline uint32_t Ptr2Indx(void *p)
	{
		return ((uint8_t *)p) - HeapNull;
	}
	inline void *Indx2Ptr(uint32_t indx)
	{
		return (void *)(HeapNull + indx);
	}
#endif // defined(_32_NORMAL) || defined(_64_EXOTIC)

	// #pragma pack(1)
	static __thread struct BLK_NODE {
		uint32_t Stamp;
		uint32_t NextIndx;
		BLK_NODE *getNext(void) const
		{
			return (BLK_NODE *)Indx2Ptr(NextIndx);
		}
		void setNext(BLK_NODE *p)
		{
			NextIndx = Ptr2Indx(p);
		}
		bool avail(void) const
		{
			return (NextIndx != 0);
		}
		void link(BLK_NODE *p)
		{
			p->NextIndx = NextIndx;
			setNext(p);
		}
		void unlink(void)
		{
			NextIndx = getNext()->NextIndx;
		}
		inline void *remove(void);
		inline void insert(void *pv, int NU);
	} __attribute__ ((packed)) BList[N_INDICES + 1];
	struct MEM_BLK : public BLK_NODE {
		uint32_t NU;
	} __attribute__ ((packed));
	// #pragma pack()

	inline void *BLK_NODE::remove(void)
	{
		BLK_NODE *p = getNext();
		unlink();
		Stamp--;
		return p;
	}
	inline void BLK_NODE::insert(void *pv, int NU)
	{
		MEM_BLK *p = (MEM_BLK *)pv;
		link(p);
		p->Stamp = ~uint32_t(0);
		p->NU = NU;
		Stamp++;
	}
	inline uint32_t U2B(uint32_t NU)
	{
		return 8 * NU + 4 * NU;
	}
	inline void SplitBlock(void *pv, uint32_t OldIndx,
						   uint32_t NewIndx)
	{
		uint32_t i, k, UDiff = Indx2Units[OldIndx] -
			Indx2Units[NewIndx];
		uint8_t *p = ((uint8_t *)pv) + U2B(Indx2Units[NewIndx]);
		if(Indx2Units[i = Units2Indx[UDiff - 1]] != UDiff) {
			k = Indx2Units[--i];
			BList[i].insert(p, k);
			p += U2B(k);
			UDiff -= k;
		}
		BList[Units2Indx[UDiff - 1]].insert(p, UDiff);
	}
	uint32_t GetUsedMemory(void)
	{
		uint32_t i, RetVal = SubAllocatorSize - (HiUnit - LoUnit) -
			(UnitsStart - pText);
		for(i = 0; i < N_INDICES; i++)
			RetVal -= U2B(Indx2Units[i] * BList[i].Stamp);
		return RetVal;
	}
	void StopSubAllocator(void)
	{
		if(SubAllocatorSize) {
			SubAllocatorSize = 0;
			delete [] HeapStart;
		}
	}
	bool StartSubAllocator(uint32_t SASize)
	{
		uint32_t t = SASize << 20U;
		if(SubAllocatorSize == t)
			return true;
		StopSubAllocator();
		if((HeapStart = new uint8_t[t]) == NULL)
			return false;
		SubAllocatorSize = t;
		return true;
	}
	inline void InitSubAllocator(void)
	{
		memset(BList, 0, sizeof(BList));
		HiUnit = (pText = HeapStart) + SubAllocatorSize;
		uint32_t Diff = U2B(SubAllocatorSize / 8 / UNIT_SIZE * 7);
		LoUnit = UnitsStart = HiUnit - Diff;
		GlueCount = GlueCount1 = 0;

#if !defined(_32_NORMAL) && !defined(_64_EXOTIC)
		HeapNull = HeapStart - 1;
#endif // !defined(_32_NORMAL) && !defined(_64_EXOTIC)
	}
	inline void GlueFreeBlocks(void)
	{
		uint32_t i, k, sz;
		MEM_BLK s0, *p, *p0, *p1;
		if(LoUnit != HiUnit)
			*LoUnit = 0;
		for((p0 = &s0)->NextIndx = i = 0; i <= N_INDICES; i++)
			while(BList[i].avail()) {
				p = (MEM_BLK *)BList[i].remove();
				if(!p->NU)
					continue;
				while((p1 = p + p->NU)->Stamp == ~uint32_t(0)) {
					p->NU += p1->NU;
					p1->NU = 0;
				}
				p0->link(p);
				p0 = p;
			}
		while(s0.avail()) {
			p = (MEM_BLK *)s0.remove();
			sz = p->NU;
			if(!sz)
				continue;
			for(; sz > 128; sz -= 128, p += 128)
				BList[N_INDICES - 1].insert(p, 128);
			if(Indx2Units[i = Units2Indx[sz - 1]] != sz) {
				k = sz - Indx2Units[--i];
				BList[k - 1].insert(p + (sz - k), k);
			}
			BList[i].insert(p, Indx2Units[i]);
		}
		GlueCount = 1 << (13 + GlueCount1++);
	}
	static void *AllocUnitsRare(uint32_t indx)
	{
		int32_t i = indx;
		do {
			if(++i == N_INDICES) {
				if(!GlueCount--) {
					GlueFreeBlocks();
					if(BList[i = indx].avail())
						return BList[i].remove();
				}
				else {
					i = U2B(Indx2Units[indx]);
					return (UnitsStart - pText > i) ?
						(UnitsStart -= i) : (NULL);
				}
			}
		} while(!BList[i].avail());
		void *RetVal = BList[i].remove();
		SplitBlock(RetVal, i, indx);
		return RetVal;
	}
	inline void *AllocUnits(uint32_t NU)
	{
		uint32_t indx = Units2Indx[NU - 1];
		if(BList[indx].avail())
			return BList[indx].remove();
		void *RetVal = LoUnit;
		LoUnit += U2B(Indx2Units[indx]);
		if(LoUnit <= HiUnit)
			return RetVal;
		LoUnit -= U2B(Indx2Units[indx]);
		return AllocUnitsRare(indx);
	}
	inline void *AllocContext(void)
	{
		if(HiUnit != LoUnit)
			return (HiUnit -= UNIT_SIZE);
		return (BList->avail()) ? (BList->remove()) :
			(AllocUnitsRare(0));
	}
	inline void UnitsCpy(void *Dest, void *Src, uint32_t NU)
	{
#if defined(_32_NORMAL) || defined(_64_NORMAL)
		uint32_t *p1 = (uint32_t *)Dest, *p2 = (uint32_t *)Src;
		do {
			p1[0] = p2[0];
			p1[1] = p2[1];
			p1[2] = p2[2];
			p1 += 3;
			p2 += 3;
		} while(--NU);
#else
		MEM_BLK *p1 = (MEM_BLK *)Dest, *p2 = (MEM_BLK *)Src;
		do {
			*p1++ = *p2++;
		} while(--NU);
#endif // defined(_32_NORMAL) || defined(_64_NORMAL)
	}
	inline void *ExpandUnits(void *OldPtr, uint32_t OldNU)
	{
		uint32_t i0 = Units2Indx[OldNU - 1],
			i1 = Units2Indx[OldNU - 1 + 1];
		if(i0 == i1)
			return OldPtr;
		void *ptr = AllocUnits(OldNU + 1);
		if(ptr) {
			UnitsCpy(ptr, OldPtr, OldNU);
			BList[i0].insert(OldPtr, OldNU);
		}
		return ptr;
	}
	inline void *ShrinkUnits(void *OldPtr, uint32_t OldNU,
							 uint32_t NewNU)
	{
		uint32_t i0 = Units2Indx[OldNU - 1],
			i1 = Units2Indx[NewNU - 1];
		if(i0 == i1)
			return OldPtr;
		if(BList[i1].avail()) {
			void *ptr = BList[i1].remove();
			UnitsCpy(ptr, OldPtr, NewNU);
			BList[i0].insert(OldPtr, Indx2Units[i0]);
			return ptr;
		}
		else {
			SplitBlock(OldPtr, i0, i1);
			return OldPtr;
		}
	}
	inline void FreeUnits(void *ptr, uint32_t NU)
	{
		uint32_t indx = Units2Indx[NU - 1];
		BList[indx].insert(ptr, Indx2Units[indx]);
	}
	inline void FreeUnit(void *ptr)
	{
		BList[((uint8_t *)ptr > UnitsStart + 128 * 1024) ?
			  (0) : (N_INDICES)].insert(ptr, 1);
	}
	inline void *MoveUnitsUp(void *OldPtr, uint32_t NU)
	{
		uint32_t indx;
		PrefetchData(OldPtr);
		if((uint8_t *)OldPtr > UnitsStart + 128 * 1024 ||
		   (BLK_NODE *)OldPtr >
		   BList[indx = Units2Indx[NU - 1]].getNext())
			return OldPtr;
		void *ptr = BList[indx].remove();
		UnitsCpy(ptr, OldPtr, NU);
		BList[N_INDICES].insert(OldPtr, Indx2Units[indx]);
		return ptr;
	}
	inline void PrepareTextArea(void)
	{
		AuxUnit = (uint8_t *)AllocContext();
		if(!AuxUnit)
			AuxUnit = UnitsStart;
		else if(AuxUnit == UnitsStart)
			AuxUnit = (UnitsStart += UNIT_SIZE);
	}
	static void ExpandTextArea(void)
	{
		BLK_NODE *p;
		uint32_t Count[N_INDICES], i = 0;
		memset(Count, 0, sizeof(Count));
		if(AuxUnit != UnitsStart) {
			if(*(uint32_t *)AuxUnit != ~uint32_t(0))
				UnitsStart += UNIT_SIZE;
			else
				BList->insert(AuxUnit, 1);
		}
		while((p = (BLK_NODE *)UnitsStart)->Stamp == ~uint32_t(0)) {
			MEM_BLK *pm = (MEM_BLK *)p;
			UnitsStart = (uint8_t *)(pm + pm->NU);
			Count[Units2Indx[pm->NU - 1]]++;
			i++;
			pm->Stamp = 0;
		}
		if(!i)
			return;
		for(p = BList + N_INDICES; p->NextIndx; p = p->getNext()) {
			while(p->NextIndx && !p->getNext()->Stamp) {
				Count[Units2Indx[((MEM_BLK *)p->getNext())->NU -
								 1]]--;
				p->unlink();
				BList[N_INDICES].Stamp--;
			}
			if(!p->NextIndx)
				break;
		}
		for(i = 0; i < N_INDICES; i++)
			for(p = BList + i; Count[i] != 0; p = p->getNext())
				while(!p->getNext()->Stamp) {
					p->unlink();
					BList[i].Stamp--;
					if(!--Count[i])
						break;
				}
	}
}
