#ifndef AC_CFT_INCLUDE_TRIVIALHEAP_HPP
#define AC_CFT_INCLUDE_TRIVIALHEAP_HPP

#include <cinttypes>

#include <array>

template <typename T, std::size_t Nm>
class TrivialHeap {
public:
    inline auto size() const { return sz; }

    inline void pop_back() {
        assert(sz > 0);
        --sz;
    }

    inline auto back() const {
        assert(sz > 0);
        return buf[sz - 1];
    }

    inline auto& back() {
        assert(sz > 0);
        return buf[sz - 1];
    }

    inline auto begin() const { return std::addressof(buf[0]); }
    inline auto begin() { return std::addressof(buf[0]); }
    inline auto end() const { return std::addressof(buf[sz]); }
    inline auto end() { return std::addressof(buf[sz]); }

    template <typename Comp>
    inline void insert(T elem, Comp comp) {
        assert(sz < Nm);
        auto* ptr = end();
        while (ptr != begin() && comp(*(ptr - 1), elem)) {
            *(ptr) = *(ptr - 1);
            --ptr;
        }
        *ptr = elem;
        ++sz;
    }

private:
    T buf[Nm];
    uint32_t sz = 0U;
};

#endif