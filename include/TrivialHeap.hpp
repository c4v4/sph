#ifndef AC_CFT_INCLUDE_TRIVIALHEAP_HPP
#define AC_CFT_INCLUDE_TRIVIALHEAP_HPP

#include <array>
#include <cinttypes>

template <typename T, std::size_t Nm, typename Comp>
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

    inline bool try_insert(T elem) {

        if (sz >= Nm) {
            if (!Comp()(back(), elem)) { return false; }
            pop_back();
        }

        auto* ptr = end();
        while (ptr != begin() && Comp()(*(ptr - 1), elem)) {
            *(ptr) = *(ptr - 1);
            --ptr;
        }
        *ptr = elem;
        ++sz;

        return true;
    }

    inline void clear() { sz = 0U; }

private:
    T buf[Nm];
    uint32_t sz = 0U;
};

#endif