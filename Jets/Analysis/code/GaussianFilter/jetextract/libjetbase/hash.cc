#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H

#include <jetbase/hash.h>

namespace jet {

	/////////////////////////////////////////////////////////////////

	// A. Appleby's MurmurHash 2.0

	// FIXME: MurmurHash can be implemented using SSE2

	uint64_t hash_murmur2_64(const void *key, const size_t len,
							 const uint64_t seed)
	{
		static const uint64_t m = 0xc6a4a7935bd1e995ULL;
		static const int r = 47;
		uint64_t h = seed ^ (len * m);
		const uint64_t *p_64 =
			reinterpret_cast<const uint64_t *>(key);
		const uint64_t *end = p_64 + (len >> 3);

		for(; p_64 != end; p_64++) {
			uint64_t k = *p_64;

			k *= m;
			k ^= k >> r;
			k *= m;
		
			h ^= k;
			h *= m;
		}

		const uint8_t *p_8 = reinterpret_cast<const uint8_t *>(p_64);

		switch(len & 7) {
		case 7:	h ^= static_cast<uint64_t>(p_8[6]) << 48;
		case 6:	h ^= static_cast<uint64_t>(p_8[5]) << 40;
		case 5:	h ^= static_cast<uint64_t>(p_8[4]) << 32;
		case 4:	h ^= static_cast<uint64_t>(p_8[3]) << 24;
		case 3:	h ^= static_cast<uint64_t>(p_8[2]) << 16;
		case 2:	h ^= static_cast<uint64_t>(p_8[1]) << 8;
		case 1:	h ^= static_cast<uint64_t>(p_8[0]);
			h *= m;
		};

		h ^= h >> r;
		h *= m;
		h ^= h >> r;

		return h;
	}

}
