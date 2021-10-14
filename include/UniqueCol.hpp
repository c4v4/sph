#ifndef SPH_INCLUDE_ROUTE_HPP_
#define SPH_INCLUDE_ROUTE_HPP_

#include <algorithm>

#include "IndexList.hpp"
#include "cft.hpp"

// https://softwareengineering.stackexchange.com/questions/402542/where-do-magic-hashing-constants-like-0x9e3779b9-and-0x9e3779b1-come-from
#define INVERSE_GOLDEN_RATIO 0x9e3779b97f4a7c15
#define LARGE_PRIME_RATIO 0x9e3779b97f4a7c55

class UniqueCol : public Column {
public:
    template <unsigned long prime>
    static inline size_t dist_hash_op(size_t seed, idx_t idx) {
        return seed ^ (idx + prime + (seed << 6) + (seed >> 2));
    }

    template <unsigned long prime>
    static inline size_t comb_hash_op(size_t seed, idx_t idx) {
        return seed ^ (idx * prime);
    }

    template <unsigned long prime>
    struct dispHash {
        size_t operator()(const Column& col) const {
            size_t seed = col.size() * col.get_cost();
            for (idx_t idx : col) { seed = dist_hash_op<prime>(seed, idx); }
            return seed;
        }
    };

    template <unsigned long prime>
    struct combHash {
        size_t operator()(const Column& col) const {
            size_t seed = col.size();
            for (idx_t idx : col) { seed = comb_hash_op<prime>(seed, idx); }
            return seed;
        }
    };

    using DispHash1 = dispHash<INVERSE_GOLDEN_RATIO>;
    using DispHash2 = dispHash<LARGE_PRIME_RATIO>;
    using CombHash1 = combHash<INVERSE_GOLDEN_RATIO>;
    using CombHash2 = combHash<LARGE_PRIME_RATIO>;

public:
    template <typename... Args>
    UniqueCol(Args&&... _args)
        : Column(std::forward<Args>(_args)...),
          disp_hash1(DispHash1()(*this)),
          disp_hash2(DispHash2()(*this)),
          comb_hash1(CombHash1()(*this)),
          comb_hash2(CombHash2()(*this)) { }

    bool operator<(const UniqueCol& other) {
        return get_cost() < other.get_cost() && size() == other.size() && comb_hash1 == other.comb_hash1 && comb_hash2 == other.comb_hash2;
    }

    bool operator>(const UniqueCol& other) {
        return get_cost() > other.get_cost() && size() == other.size() && comb_hash1 == other.comb_hash1 && comb_hash2 == other.comb_hash2;
    }

    bool operator==(const UniqueCol& other) {
        return get_cost() == other.get_cost() && size() == other.size() && disp_hash1 == other.disp_hash1 && disp_hash2 == other.disp_hash2 &&
               comb_hash1 == other.comb_hash1 && comb_hash2 == other.comb_hash2;
    }

    [[nodiscard]] inline size_t get_disp_hash1() const { return disp_hash1; }
    [[nodiscard]] inline size_t get_disp_hash2() const { return disp_hash2; }
    [[nodiscard]] inline size_t get_comb_hash1() const { return comb_hash1; }
    [[nodiscard]] inline size_t get_comb_hash2() const { return comb_hash2; }

private:
    size_t disp_hash1;
    size_t disp_hash2;
    size_t comb_hash1;
    size_t comb_hash2;
};

#endif