#ifndef SPH_INCLUDE_ROUTEPOOL_HPP_
#define SPH_INCLUDE_ROUTEPOOL_HPP_

#include <cassert>

#include "UniqueCol.hpp"
#include "VectorSet.hpp"

namespace sph {

    class UniqueColSet {
    private:
        struct Hash {
            size_t operator()(const UniqueCol &column) const { return column.size() ^ column.get_comb_hash1() ^ column.get_comb_hash2(); }
        };

    public:
        UniqueColSet() { }

        template <typename Iter>
        UniqueColSet(Iter beg, Iter end) : vec_set(beg, end) { }

        void reserve(idx_t size) { vec_set.reserve(size); }

        template <typename... _Args>
        std::pair<idx_t, bool> add_column(_Args &&...args) {

            UniqueCol candidate_elem(std::forward<_Args>(args)...);

            auto [set_elem, inserted_new] = vec_set.emplace_back(candidate_elem);
            assert(inserted_new ||
                   (candidate_elem.size() == set_elem->size() && candidate_elem.get_comb_hash1() == set_elem->get_comb_hash1() &&
                    candidate_elem.get_comb_hash2() == set_elem->get_comb_hash2()));

            if (!inserted_new) {                                               // ==> same combination and size
                if (candidate_elem.get_solcost() < set_elem->get_solcost() ||  // Keep the one belonging to the best sol ...
                    candidate_elem < *set_elem) {                              // ... or the dominant one
                    *set_elem = std::move(candidate_elem);
                }
            }
            return std::make_pair(set_elem - vec_set.begin(), inserted_new);
        }

        template <typename Iter>
        size_t insert(Iter beg, Iter end) {
            return std::accumulate(beg, end, 0, [this](size_t sum, auto it) { return sum + add_column(UniqueCol(*it)).second; });
        }

        [[nodiscard]] inline size_t size() const { return vec_set.size(); }
        [[nodiscard]] inline auto begin() { return vec_set.begin(); }
        [[nodiscard]] inline auto end() { return vec_set.end(); }
        [[nodiscard]] inline auto begin() const { return vec_set.begin(); }
        [[nodiscard]] inline auto end() const { return vec_set.end(); }
        [[nodiscard]] inline UniqueCol &operator[](size_t i) { return vec_set[i]; }
        [[nodiscard]] inline const UniqueCol &operator[](size_t i) const { return vec_set[i]; }

    private:
        cav::VectorSet<UniqueCol, Hash> vec_set;
    };

}  // namespace sph

#endif