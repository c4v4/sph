#ifndef SPH_INCLUDE_INDEXLIST_HPP_
#define SPH_INCLUDE_INDEXLIST_HPP_

#include <algorithm>
#include <vector>

#include "cft.hpp"

class IndexList {
public:
    IndexList() : sz(0U), data(nullptr) { }
    IndexList(const idx_t* beg_, const idx_t* end_) : sz(end_ - beg_), data(new idx_t[sz]) { std::copy(beg_, end_, data); }
    IndexList(const IndexList& other) : IndexList(other.begin(), other.end()) { }
    IndexList(IndexList&& other) : sz(other.sz), data(std::exchange(other.data, nullptr)) { }
    IndexList(idx_t*& other_data, idx_t other_sz) : sz(other_sz), data(std::exchange(other_data, nullptr)) { }

    IndexList& operator=(const IndexList& other) {
        if (sz < other.sz) {
            delete[] data;
            data = new idx_t[other.sz];
        }
        sz = other.sz;
        std::copy(other.begin(), other.end(), data);
        return *this;
    };

    IndexList& operator=(IndexList&& other) {
        delete[] data;
        data = std::exchange(other.data, nullptr);
        sz = std::exchange(other.sz, 0);
        return *this;
    };

    ~IndexList() { delete[] data; };

    [[nodiscard]] inline idx_t* begin() { return data; }
    [[nodiscard]] inline idx_t* end() { return begin() + sz; }
    [[nodiscard]] inline idx_t& operator[](idx_t i) { return data[i]; }
    [[nodiscard]] inline idx_t& front() { return data[0]; }
    [[nodiscard]] inline idx_t& back() { return data[sz - 1]; }

    [[nodiscard]] inline idx_t size() const { return sz; }
    [[nodiscard]] inline bool empty() const { return sz == 0; }
    [[nodiscard]] inline const idx_t* begin() const { return data; }
    [[nodiscard]] inline const idx_t* end() const { return begin() + sz; }
    [[nodiscard]] inline idx_t operator[](idx_t i) const { return data[i]; }
    [[nodiscard]] inline idx_t front() const { return data[0]; }
    [[nodiscard]] inline idx_t back() const { return data[sz - 1]; }

private:
    idx_t sz;
    idx_t* data;
};


class Column : public IndexList {

public:
    template <typename Col>
    Column(const Col& col) : Column(col.begin(), col.end(), col.get_cost(), col.get_solcost()) { }

    template <typename Iter>
    Column(Iter beg_, Iter end_, real_t c_, real_t sol_c_) : IndexList(beg_, end_), c(c_), sol_c(sol_c_) { }

    Column(idx_t*& other_data, idx_t other_sz, real_t c_, real_t sol_c_) : IndexList(other_data, other_sz), c(c_), sol_c(sol_c_) { }
    Column(const Column& other) : IndexList(other.begin(), other.end()), c(other.c), sol_c(other.sol_c) { }
    Column(Column&& other) : IndexList(std::move(other)), c(std::exchange(other.c, 0)), sol_c(std::exchange(other.sol_c, 0)) { }

    Column& operator=(const Column& other) {
        IndexList::operator=(other);
        c = other.c;
        sol_c = other.sol_c;
        return *this;
    };

    Column& operator=(Column&& other) {
        IndexList::operator=(std::move(other));
        c = std::exchange(other.c, 0);
        sol_c = std::exchange(other.sol_c, 0);
        return *this;
    };

    [[nodiscard]] inline auto get_cost() const { return c; }

    [[nodiscard]] inline auto get_solcost() const { return sol_c; }
    inline void set_solcost(real_t new_c) { sol_c = new_c; }


    template <typename Multipliers>
    [[nodiscard]] inline auto compute_lagr_cost(const Multipliers& u) const {
        real_t local_c_u = c;
        for (idx_t i : *this) {
            assert(i < u.size());
            local_c_u -= u[i];
        }
        return local_c_u;
    }

    bool operator==(const Column& other) const {
        if (c != other.c || sol_c != other.sol_c || size() != other.size()) { return false; }
        for (idx_t n = 0; n < size(); ++n) {
            if ((*this)[n] != other[n]) { return false; }
        }
        return true;
    }

private:
    real_t c = 0.0;
    real_t sol_c = 0.0;
};


class Row : public std::vector<idx_t> {
public:
    using std::vector<idx_t>::vector;
};


#endif