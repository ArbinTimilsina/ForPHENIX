namespace ppmd {

	enum { TOP = 1 << 24, BOT = 1 << 15 };
	static __thread struct SUBRANGE {
		uint32_t low, high, scale;
	} Range;
	static __thread uint32_t low, code, range;

	inline void rcInitEncoder(void)
	{
		low = 0;
		range = uint32_t(-1);
	}
	template<typename Stream_t>
	inline void RC_ENC_NORMALIZE(Stream_t stream)
	{
		while((low ^ (low + range)) < TOP || range < BOT &&
			  ((range = -low & (BOT - 1)), 1)) {
			_PPMD_E_PUTC(low >> 24, stream);
			range <<= 8;
			low <<= 8;
		}
	}
	inline void rcEncodeSymbol(void)
	{
		low += Range.low * (range /= Range.scale);
		range *= Range.high - Range.low;
	}
	template<typename Stream_t>
	inline void rcFlushEncoder(Stream_t stream)
	{
		for(uint32_t i = 0; i < 4; i++) {
			_PPMD_E_PUTC(low >> 24, stream);
			low <<= 8;
		}
	}
	template<typename Stream_t>
	inline void rcInitDecoder(Stream_t stream)
	{
		low = code = 0;
		range = uint32_t(-1);
		for(uint32_t i = 0; i < 4; i++)
			code = (code << 8) | _PPMD_D_GETC(stream);
	}
	template<typename Stream_t>
	inline void RC_DEC_NORMALIZE(Stream_t stream)
	{
		while((low ^ (low + range)) < TOP || range < BOT &&
			  ((range = -low & (BOT - 1)), 1)) {
			code = (code << 8) | _PPMD_D_GETC(stream);
			range <<= 8;
			low <<= 8;
		}
	}
	inline uint32_t rcGetCurrentCount(void)
	{
		return (code - low) / (range /= Range.scale);
	}
	inline void rcRemoveSubrange(void)
	{
		low += range * Range.low;
		range *= Range.high - Range.low;
	}
	inline uint32_t rcBinStart(uint32_t f0, uint32_t Shift)
	{
		return f0 * (range >>= Shift);
	}
	inline uint32_t rcBinDecode(uint32_t tmp)
	{
		return (code - low >= tmp);
	}
	inline void rcBinCorrect0(uint32_t tmp)
	{
		range = tmp;
	}
	inline void rcBinCorrect1(uint32_t tmp, uint32_t f1)
	{
		low += tmp;
		range *= f1;
	}
}
