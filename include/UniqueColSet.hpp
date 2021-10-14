#ifndef SPH_INCLUDE_ROUTEPOOL_HPP_
#define SPH_INCLUDE_ROUTEPOOL_HPP_

#include <cassert>

#include "UniqueCol.hpp"
#include "VectorSet.hpp"

class UniqueColSet {
private:
    struct Hash {
        size_t operator()(const UniqueCol &column) const {
            return column.get_disp_hash1() ^ column.get_disp_hash2() ^ column.get_comb_hash1() ^ column.get_comb_hash2();
        }
    };

public:
    UniqueColSet() { }

    template <typename Iter>
    UniqueColSet(Iter beg, Iter end) : vec_set(beg, end) { }

    void reserve(idx_t size) { vec_set.reserve(size); }

    template <typename... _Args>
    bool add_column(_Args &&...args) {
        UniqueCol candidate_elem(std::forward<_Args>(args)...);
        auto [set_elem, inserted] = vec_set.emplace_back(candidate_elem);

        if (!inserted && candidate_elem < *set_elem) {
            *set_elem = std::move(candidate_elem);
            return false;
        }

        return true;
    }

    template <typename Iter>
    size_t insert(Iter beg, Iter end) {
        return std::accumulate(beg, end, 0, [this](size_t sum, auto it) { return sum + add_column(UniqueCol(*it)); });
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

#endif