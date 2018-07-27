#ifndef XJETSTAT_MIXING_H_
#define XJETSTAT_MIXING_H_

#include <jetstat/histogram.h>
#include <deque>

namespace jet {

	class event_embedding_pool_t {
	private:
		static const unsigned int nbin_centrality = 5;
		static const unsigned int nbin_vertex = 10;
		class life_time_t {
		protected:
			uint64_t _usage;
			event_t _content;
		public:
			inline life_time_t(
				const uint64_t usage, const event_t content)
				: _usage(usage), _content(content)
			{
			}
			inline const event_t &operator()(void)
			{
				_usage--;

				return _content;
			}
			inline bool expired(void) const
			{
				return _usage <= 0;
			}
		};
		bool _second_is_empty;
		std::deque<life_time_t> _first[(nbin_centrality + 1) * nbin_vertex];
		std::deque<life_time_t> _second[(nbin_centrality + 1) * nbin_vertex];
	private:
		inline void prune(std::deque<life_time_t> &deque)
		{
			deque.erase(
				std::remove_if(
					deque.begin(), deque.end(),
					std::mem_fun_ref(&life_time_t::expired)),
				deque.end());
		}
		inline void prune(std::deque<life_time_t> *histogram)
		{
			for(int i = 0; i < (nbin_centrality + 1) * nbin_vertex; i++) {
				prune(_first[i]);
				prune(_second[i]);
			}
		}
		inline void prune(void)
		{
			prune(_first);
			prune(_second);
		}
		inline int flat_index(const event_t &event)
		{
			if(!(event.vertex_height() >= -25.0F &&
				 event.vertex_height() < 25.0F)) {
				return -1;
			}

			const int index_vertex_height = static_cast<int>(
				floorf((event.vertex_height() + 25.0F) *
					   (nbin_vertex / 50.0F)));
			const int index_centrality =
				(event.centrality() >= 1 &&
				 event.centrality() < 101) ?
				(static_cast<int>(rintf(
					event.centrality())) - 1) / 20 :
				nbin_centrality;
			return index_vertex_height * (nbin_centrality + 1) +
				index_centrality;
		}
	public:
		inline event_embedding_pool_t(const bool second_is_empty)
			: _second_is_empty(second_is_empty)
		{
		}
		inline void fill_first(
			const uint64_t usage, const event_t &entry)
		{
			int index = flat_index(entry);

			if(index >= 0) {
				_first[index].push_back(
					life_time_t(usage, entry));
			}
		}
		inline void fill_second(
			const uint64_t usage, const event_t &entry)
		{
			int index = flat_index(entry);

			if(index >= 0) {
				_second[index].push_back(
					life_time_t(usage, entry));
			}
		}
		inline bool pop_first(event_t &entry)
		{
			for(int i = 0; i < (nbin_centrality + 1) * nbin_vertex; i++) {
				if(!_first[i].empty()) {
					entry = (_first[i].front())();
					prune();
					return true;
				}
			}
			return false;
		}
		inline bool pop_pair(std::pair<event_t, event_t> &entry)
		{
			if(_second_is_empty) {
				return pop_first(entry.first);
			}

			for(int i = 0; i < (nbin_centrality + 1) * nbin_vertex; i++) {
				if(!_first[i].empty() && !_second[i].empty()) {
					entry.first = (_first[i].front())();
					entry.second = (_second[i].front())();
					prune();
					return true;
				}
			}
			return false;
		}
		inline size_t embeddable(void) const
		{
			for(int i = 0; i < (nbin_centrality + 1) * nbin_vertex; i++) {
				if(!_first[i].empty() && _second[i].empty()) {
					return false;
				}
			}

			return true;
		}
	};

}

#endif // XJETSTAT_MIXING_H_
