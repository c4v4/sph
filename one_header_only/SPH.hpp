/* Automatically generated one-header-only library */


/* ################################################################# */
/* #### Original Header: include/cft.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_CFT_HPP_
#define SPH_INCLUDE_CFT_HPP_

#include <cstdint>
#include <limits>

/* Types and constants */
typedef uint32_t idx_t;
typedef double real_t;

#define REAL_MAX (std::numeric_limits<real_t>::max())
#define REAL_LOWEST (std::numeric_limits<real_t>::lowest())
#define REMOVED_INDEX (std::numeric_limits<idx_t>::max())
#define MAX_INDEX (std::numeric_limits<idx_t>::max())
#define HAS_INTEGRAL_COSTS 1.0  // 1.0 if yes , 0.0 if no

#ifndef INST_HARD_CAP
#define INST_HARD_CAP 200'000U
#endif

/* UTILS */
#define NO_INLINE  //__attribute__((noinline))

#ifdef VERBOSE
#define IF_VERBOSE
#else
#define IF_VERBOSE if constexpr (false)
#endif

#ifdef NDEBUG
#define IF_DEBUG if constexpr (false)
#else
#define IF_DEBUG
#endif

#endif  // SPH_INCLUDE_CFT_HPP_


/* ################################################################# */
/* #### Original Header: include/CollectionOf.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_COLUMNS_HPP_
#define SPH_INCLUDE_COLUMNS_HPP_
#include <fmt/core.h>
#include <stddef.h>

#include <algorithm>
#include <cassert>
#include <numeric>

/* #include "cft.hpp" */

#define GROWTH_FACTOR 1.618

template <typename Elem>
class CollectionOf;

template <typename Elem, class Ptr>
static Ptr* align_ptr(Ptr* src) {

    size_t alignment = reinterpret_cast<uintptr_t>(reinterpret_cast<const void*>(src)) % alignof(Elem);
    if (alignment > 0) { src += alignof(Elem) - alignment; }
    return src;
}

template <typename Elem>
class IdxList {
public:
    static constexpr size_t MY_SIZE = sizeof(Elem);

protected:
    IdxList(const idx_t* beg_, const idx_t* end_) : sz(end_ - beg_) { std::copy(beg_, end_, begin()); }
    IdxList(const IdxList& other) : sz(other.sz) { std::copy(other.begin(), other.end(), begin()); }

public:
    IdxList() : sz(0U) { }
    IdxList(IdxList&&) = delete;
    IdxList& operator=(const IdxList&) = delete;
    IdxList& operator=(IdxList&&) = delete;
    ~IdxList() = default;

    Elem& operator=(const Elem&) = delete;
    Elem& operator=(Elem&&) = delete;

    [[nodiscard]] inline idx_t* begin() { return reinterpret_cast<idx_t*>(reinterpret_cast<char*>(this) + MY_SIZE); }
    [[nodiscard]] inline idx_t* end() { return begin() + sz; }
    [[nodiscard]] inline idx_t& operator[](idx_t i) { return begin()[i]; }

    [[nodiscard]] inline idx_t size() const { return sz; }
    [[nodiscard]] inline bool empty() const { return sz == 0; }
    [[nodiscard]] inline const idx_t* begin() const { return reinterpret_cast<const idx_t*>(reinterpret_cast<const char*>(this) + MY_SIZE); }
    [[nodiscard]] inline const idx_t* end() const { return begin() + sz; }
    [[nodiscard]] inline idx_t operator[](idx_t i) const { return begin()[i]; }

protected:
    idx_t sz;
};

class SubInstCol : public IdxList<SubInstCol> {
    friend class CollectionOf<SubInstCol>;

protected:
    SubInstCol(const idx_t* beg_, const idx_t* end_, real_t c_ = 0.0) : IdxList<SubInstCol>(beg_, end_), c(c_), c_u(c_) { }
    SubInstCol(const SubInstCol& other) : IdxList<SubInstCol>(other), c(other.c), c_u(other.c_u) { }

public:
    SubInstCol(real_t c_ = 0.0) : c(c_), c_u(c_) { }

    [[nodiscard]] inline real_t get_cost() const { return c; }
    inline void set_cost(real_t new_c) { c_u = c = new_c; }

    [[nodiscard]] inline real_t get_cu() const { return c_u; }
    inline void set_cu(real_t new_c) { c_u = new_c; }
    inline real_t update_cu(const real_t delta_u) { return c_u -= delta_u; }


    [[nodiscard]] inline real_t compute_lagr_cost(const std::vector<real_t>& u) const {
        real_t local_c_u = c;
        for (idx_t i : *this) {
            assert(i < u.size());
            local_c_u -= u[i];
        }
        return local_c_u;
    }

    inline real_t compute_lagr_cost(const std::vector<real_t>& u) {
        c_u = c;
        for (idx_t i : *this) {
            assert(i < u.size());
            c_u -= u[i];
        }
        return c_u;
    }

protected:
    real_t c;
    real_t c_u;
};

static_assert(sizeof(SubInstCol) == SubInstCol::MY_SIZE);
static_assert(SubInstCol::MY_SIZE == IdxList<SubInstCol>::MY_SIZE);
static_assert(IdxList<SubInstCol>::MY_SIZE == std::max(alignof(real_t), sizeof(idx_t)) + sizeof(real_t) * 2U);


template <typename Elem>
class CollectionOf {
public:
    class CollectionIter {
    public:
        CollectionIter(Elem* base_) : base(reinterpret_cast<Elem*>(align_ptr<Elem>(reinterpret_cast<char*>(base_)))) {
            assert(reinterpret_cast<uintptr_t>(reinterpret_cast<const void*>(base)) % alignof(Elem) == 0);
        }

        [[nodiscard]] inline auto& operator*() { return *base; }

        inline auto& operator++() {
            base = reinterpret_cast<Elem*>(align_ptr<Elem>(reinterpret_cast<char*>(base) + sizeof(Elem) + base->size() * sizeof(idx_t)));
            assert(align_ptr<Elem>(base) == base);
            return *this;
        }

        inline auto operator++(int) {
            Elem* old_base = base;
            base = reinterpret_cast<Elem*>(align_ptr<Elem>(reinterpret_cast<char*>(base) + sizeof(Elem) + base->size() * sizeof(idx_t)));
            assert(align_ptr<Elem>(base) == base);
            return CollectionIter(old_base);
        }

        [[nodiscard]] inline bool operator!=(CollectionIter& it2) const { return base != it2.base; }

        [[nodiscard]] inline bool operator==(CollectionIter& it2) const { return base == it2.base; }

        [[nodiscard]] Elem* data() const { return base; }

    private:
        Elem* base;
    };


public:
    CollectionOf(idx_t allocated_elem = 100) {
        buffer_allocate(allocated_elem * (sizeof(Elem) + sizeof(idx_t) * 10));
        offsets.reserve(allocated_elem);
    }

    ~CollectionOf() {
        free(start);
        start = finish = end_of_storage = nullptr;
    }

    [[nodiscard]] inline CollectionIter begin() { return CollectionIter(reinterpret_cast<Elem*>(start)); }
    [[nodiscard]] inline CollectionIter end() { return CollectionIter(reinterpret_cast<Elem*>(finish)); }
    [[nodiscard]] inline Elem& operator[](idx_t j) { return *reinterpret_cast<Elem*>(start + offsets[j]); }
    [[nodiscard]] inline Elem& back() { return *reinterpret_cast<Elem*>(start + offsets.back()); }
    [[nodiscard]] inline Elem* data() { return reinterpret_cast<Elem*>(start + offsets[0]); }

    [[nodiscard]] inline idx_t size() const { return offsets.size(); }
    [[nodiscard]] inline bool empty() const { return offsets.empty(); }
    [[nodiscard]] inline const CollectionIter begin() const { return CollectionIter(reinterpret_cast<Elem*>(start)); }
    [[nodiscard]] inline const CollectionIter end() const { return CollectionIter(reinterpret_cast<Elem*>(finish)); }
    [[nodiscard]] inline const Elem& operator[](idx_t j) const { return *reinterpret_cast<Elem*>(start + offsets[j]); }
    [[nodiscard]] inline const Elem& back() const { return *reinterpret_cast<Elem*>(start + offsets.back()); }
    [[nodiscard]] inline const Elem* data() const { return reinterpret_cast<Elem*>(align_ptr<Elem>(start)); }

    inline void clear() {
        finish = start;
        offsets.clear();
        assert(!is_corrupted());
    }

    inline void reserve(idx_t new_ncols) {
        idx_t ncols = offsets.size();
        if (new_ncols <= ncols) { return; }

        double avg_col_size = ncols > 0 ? static_cast<double>(finish - start) / static_cast<double>(ncols) : 10.0;

        buffer_allocate(new_ncols * avg_col_size + 1);
        offsets.reserve(new_ncols);
        assert(!is_corrupted());
    }

    template <typename iter, typename... Args>
    inline Elem& emplace_back(iter beg_, iter end_, Args&&... args) {

        idx_t col_size = end_ - beg_;
        idx_t col_storage = sizeof(Elem) + col_size * sizeof(idx_t);
        char* new_finish = align_ptr<Elem>(finish) + col_storage;

        buffer_check_size(new_finish);
        Elem* elem = new (align_ptr<Elem>(finish)) Elem(beg_, end_, std::forward<Args>(args)...);

        offsets.push_back(align_ptr<Elem>(finish) - start);
        finish = new_finish;

        assert(!is_corrupted());

        return *elem;
    }

    inline Elem& emplace_back(const Elem& other) {

        idx_t col_size = other.size();
        idx_t col_storage = sizeof(Elem) + col_size * sizeof(idx_t);
        char* new_finish = align_ptr<Elem>(finish) + col_storage;

        buffer_check_size(new_finish);
        Elem* elem = new (align_ptr<Elem>(finish)) Elem(other);

        offsets.push_back(align_ptr<Elem>(finish) - start);
        finish = new_finish;

        assert(!is_corrupted());

        return *elem;
    }

    inline void push_back(const Elem& elem) { emplace_back(elem); }


    // Construct in place a new column
    template <typename... Args>
    inline Elem& new_col_create(Args&&... args) {
        Elem& elem = buffer_emplace_back<Elem>(std::forward<Args>(args)...);
        offsets.push_back(reinterpret_cast<char*>(std::addressof(elem)) - start);
        assert(!is_corrupted());
        return elem;
    }

    inline void new_col_push_back(idx_t idx) {
        buffer_emplace_back<idx_t>(idx);
        ++back().sz;
        assert(!is_corrupted());
    }

    inline void new_col_discard() {
        finish = reinterpret_cast<char*>(offsets.back());
        offsets.pop_back();
        assert(!is_corrupted());
    }


private:
    template <typename T, typename... Args>
    inline T& buffer_emplace_back(Args&&... args) {
        char* new_finish = align_ptr<T>(finish) + sizeof(T);

        buffer_check_size(new_finish);
        T* elem = new (align_ptr<T>(finish)) T(std::forward<Args>(args)...);

        finish = new_finish;

        return *elem;
    }

    inline void buffer_check_size(char*& new_finish) {
        if (new_finish > end_of_storage) {
            size_t new_size = std::max<size_t>(new_finish - start, GROWTH_FACTOR * (end_of_storage - start));
            char* old_start = start;
            buffer_allocate(new_size);
            new_finish += start - old_start;
        }
    }

    inline void buffer_allocate(size_t new_size) {
        assert(!start == !finish && !start == !end_of_storage);

        size_t prev_size = finish - start;
        start = reinterpret_cast<char*>(realloc(start, new_size));
        finish = start + prev_size;
        end_of_storage = start + new_size;

        assert(start && finish && end_of_storage);
    }

#define CHECK_RET is_corrupted_printed()
    //#define CHECK_RET true

    bool is_corrupted() {

        if (align_ptr<Elem>(begin().data()) != begin().data()) { return CHECK_RET; }
        if (!empty() && align_ptr<Elem>(begin().data()) != reinterpret_cast<Elem*>(start + offsets[0])) { return CHECK_RET; }
        if (align_ptr<Elem>(end().data()) != end().data()) { return CHECK_RET; }

        idx_t counter = 0;
        for (Elem& elem : *this) {
            if (align_ptr<Elem>(std::addressof(elem)) != std::addressof(elem)) { return CHECK_RET; }
            if (++counter > size()) { return CHECK_RET; }
        }
        return false;
    }

    bool is_corrupted_printed() {

        fmt::print("begin {} {}\n", (void*)begin().data(), (void*)align_ptr<Elem>(begin().data()));
        if (align_ptr<Elem>(begin().data()) != begin().data()) { return true; }

        if (!empty()) { fmt::print("offset begin {} {}\n", (void*)begin().data(), (void*)(start + offsets[0])); }
        if (!empty() && align_ptr<Elem>(begin().data()) != reinterpret_cast<Elem*>(start + offsets[0])) { return true; }

        fmt::print("end {} {}\n", (void*)end().data(), (void*)align_ptr<Elem>(end().data()));
        if (align_ptr<Elem>(end().data()) != end().data()) { return true; }

        idx_t counter = 0;
        for (Elem& elem : *this) {
            fmt::print("{} ", (void*)std::addressof(elem));
            if (align_ptr<Elem>(std::addressof(elem)) != std::addressof(elem)) { return true; }
            if (++counter > size()) {
                fmt::print("\n");
                return true;
            }
        }
        fmt::print("\n");

        return false;
    }

    char* start = nullptr;
    char* finish = nullptr;
    char* end_of_storage = nullptr;

    std::vector<size_t> offsets;
};


// typedef CollectionOf<RowTempName> Rows;
// typedef CollectionOf<InstCol> Cols;
typedef CollectionOf<SubInstCol> SubInstCols;

#endif

/* ################################################################# */
/* #### Original Header: include/MStar.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_MSTAR_HPP_
#define SPH_INCLUDE_MSTAR_HPP_

/* #include "cft.hpp" */
class MStar {
public:
    MStar() = default;

    explicit MStar(idx_t nrows) : M_star(nrows, 0), zeros(nrows) { }

    inline void reset_uncovered(idx_t nrows) {
        zeros = nrows;
        M_star.assign(nrows, 0);
    }

    template<typename Collection>
    inline void reset_covered(const Collection& cols, idx_t nrows) {
        reset_uncovered(nrows);
        for (auto& j_col : cols) { cover_rows(j_col); }
    }

    template<typename Collection>
    inline void reset_covered(const Collection& cols, const std::vector<idx_t>& indexes, idx_t nrows) {
        reset_uncovered(nrows);
        for (auto& j : indexes) { cover_rows(cols[j]); }
    }

    inline bool is_redundant(const SubInstCol& col) {
        for (const auto i : col) {
            assert(M_star[i] > 0);
            if (M_star[i] == 1) { return false; }
        }
        return true;
    }

    template<typename IterableList>
    inline void cover_rows(const IterableList& rows) {
        for (auto i : rows) {
            assert(i < M_star.size());
            cover(i);
        }
    }

    template<typename IterableList>
    inline void uncover_rows(const IterableList& rows) {
        for (auto i : rows) {
            assert(i < M_star.size());
            uncover(i);
        }
    }

    inline void cover(idx_t row) {
        zeros -= static_cast<idx_t>(M_star[row] == 0);
        ++M_star[row];
    }

    inline void uncover(idx_t row) {
        zeros += static_cast<idx_t>(M_star[row] == 1);
        --M_star[row];
    }

    [[nodiscard]] idx_t get(idx_t idx) const { return M_star[idx]; }
    idx_t operator[](idx_t idx) const { return M_star[idx]; }
    idx_t& operator[](idx_t idx) { return M_star[idx]; }

    [[nodiscard]] inline auto begin() const { return M_star.begin(); }
    [[nodiscard]] inline auto end() const { return M_star.end(); }
    inline auto begin() { return M_star.begin(); }
    inline auto end() { return M_star.end(); }

    [[nodiscard]] idx_t get_uncovered() const { return zeros; }
    [[nodiscard]] idx_t get_covered() const { return M_star.size() - zeros; }

    [[nodiscard]] idx_t size() const { return M_star.size(); }


private:
    std::vector<idx_t> M_star;
    idx_t zeros{};
};

#endif

/* ################################################################# */
/* #### Original Header: include/Stopwatch.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_STOPWATCH_HPP_
#define SPH_INCLUDE_STOPWATCH_HPP_

#include <chrono>

class Stopwatch {
public:
    Stopwatch() : begin(std::chrono::steady_clock::now()), last(begin) { }

    double seconds_from_begin() { return std::chrono::duration<double>(now() - begin).count(); }
    double millisec_from_begin() { return std::chrono::duration<double, std::milli>(now() - begin).count(); }
    double microsec_from_begin() { return std::chrono::duration<double, std::micro>(now() - begin).count(); }
    double nanosec_from_begin() { return std::chrono::duration<double, std::nano>(now() - begin).count(); }

    double seconds_lap() {
        auto old_last = last;
        last = now();
        return std::chrono::duration<double>(last - old_last).count();
    }

    double millisec_lap() {
        auto old_last = last;
        last = now();
        return std::chrono::duration<double, std::milli>(last - old_last).count();
    }

    double microsec_lap() {
        auto old_last = last;
        last = now();
        return std::chrono::duration<double, std::micro>(last - old_last).count();
    }

    double nanosec_lap() {
        auto old_last = last;
        last = now();
        return std::chrono::duration<double, std::nano>(last - old_last).count();
    }

protected:
    std::chrono::steady_clock::time_point now() { return std::chrono::steady_clock::now(); }

private:
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point last;
};


class Timer : public Stopwatch {

public:
    Timer() : end(now()) { }
    Timer(double second_tlim) : end(now() + std::chrono::microseconds(static_cast<uint64_t>(second_tlim * 1E6))) { }

    double seconds_until_end() { return std::chrono::duration<double>(end - now()).count(); }
    double millisec_until_end() { return std::chrono::duration<double, std::milli>(end - now()).count(); }
    double microsec_until_end() { return std::chrono::duration<double, std::micro>(end - now()).count(); }
    double nanosec_until_end() { return std::chrono::duration<double, std::nano>(end - now()).count(); }

    bool exceeded_tlim() { return now() >= end; }

private:
    std::chrono::steady_clock::time_point end;
};

#endif

/* ################################################################# */
/* #### Original Header: include/IndexList.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_INDEXLIST_HPP_
#define SPH_INCLUDE_INDEXLIST_HPP_

#include <algorithm>
#include <vector>

/* #include "cft.hpp" */

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

/* ################################################################# */
/* #### Original Header: include/VectorSet.hpp */
/* ################################################################# */

#ifndef CAV_VECTORSET_HPP
#define CAV_VECTORSET_HPP

#include <type_traits>
#include <unordered_set>
#include <vector>

namespace cav {

    template <typename T, typename Hash = std::hash<T>>
    class VectorSet {
    private:
        using viter = typename std::vector<T>::iterator;
        using const_viter = typename std::vector<T>::const_iterator;

        /**
         * @brief Since the set does not contain T, but *T,
         * this hacky struct is used to make the *T behave like
         * a T would do (for what the set is concerned).
         *
         * Note: it's still a ptr, so it has all the problems of
         * a ptr.
         * NEVER return it if it is initialized to a local var.
         */
        struct TPtrWrap {
            TPtrWrap(T*& base_addr_, size_t idx_) : base_addr(base_addr_), idx(idx_) { }
            operator T&() const { return base_addr[idx]; }
            T* data() const { return base_addr + idx; }
            bool operator==(const TPtrWrap& other) const { return other.base_addr[other.idx] == base_addr[idx]; }

        private:
            T*& base_addr;
            size_t idx;
        };

        using Set = std::unordered_set<TPtrWrap, Hash>;
        using siter = typename Set::iterator;
        using const_siter = typename Set::const_iterator;

    public:
        VectorSet() { }

        template <typename Iter>
        VectorSet(Iter beg, Iter end) {
            reserve(end - beg);
            for (auto it = beg; it != end; ++it) { emplace_back(*it); }
        }

        VectorSet(const VectorSet& other) : vec(other.vec), set(other.set) { vec_data = vec.data(); }

        VectorSet(VectorSet&& other) : vec(std::move(other.vec)), set(std::move(other.set)) { vec_data = vec.data(); }

        VectorSet& operator=(const VectorSet& other) {
            vec = other.vec;
            set = other.set;
            vec_data = vec.data();
            return *this;
        }

        VectorSet& operator=(VectorSet&& other) {
            vec = std::move(other.vec);
            set = std::move(other.set);
            vec_data = vec.data();
            other.vec_data = nullptr;
            return *this;
        }

        template <typename... _Args>
        std::pair<viter, bool> emplace_back(_Args&&... args) {
            T elem(std::forward<_Args>(args)...);
            T* elem_ptr = std::addressof(elem);
            auto wrapped_elem = TPtrWrap(elem_ptr, 0);

            auto set_elem_iter = set.find(wrapped_elem);
            if (set_elem_iter == set.end()) {
                T& inserted = vec.emplace_back(std::move(elem));
                vec_data = vec.data();
                set.emplace(vec_data, vec.size() - 1);
                return std::make_pair<viter, bool>(viter(std::addressof(inserted)), true);
            }

            return std::make_pair<viter, bool>(viter(set_elem_iter->data()), false);
        }

        template <typename Tt>
        std::pair<viter, bool> push_back(Tt&& elem) {
            return emplace_back(std::forward<Tt>(elem));
        }

        void reserve(size_t size) {
            vec.reserve(size);
            set.reserve(size);
        }

        // Vector operations
        T& back() { return vec.back(); }
        T& front() { return vec.front(); }
        viter begin() { return vec.begin(); }
        viter end() { return vec.end(); }
        T& operator[](size_t i) { return vec[i]; }
        void clear() { vec.clear(), set.clear(); }

        const T& back() const { return vec.back(); }
        const T& front() const { return vec.front(); }
        const_viter begin() const { return vec.begin(); }
        const_viter end() const { return vec.end(); }
        const T& operator[](size_t i) const { return vec[i]; }
        bool empty() const { return vec.empty(); }
        size_t size() const { return vec.size(); }

        // Set Operations
        const_viter find(T elem) const {
            T* elem_ptr = std::addressof(elem);
            auto wrapped_elem = TPtrWrap(elem_ptr, 0);

            const_siter found = set.find(wrapped_elem);
            if (found != set.cend()) return const_viter(found->data());
            return vec.cend();
        }

        viter find(T elem) {
            T* elem_ptr = std::addressof(elem);
            TPtrWrap wrapped_elem(elem_ptr, 0);

            siter found = set.find(wrapped_elem);
            if (found != set.end()) return viter(found->data());
            return vec.end();
        }

        size_t count(T& elem) const { return find(elem) != vec.cend(); }
        size_t count(T& elem) { return find(elem) != vec.end(); }

        std::vector<T>& get_vector() { return vec; }

    private:
        std::vector<T> vec;
        T* vec_data = nullptr;  // we can't get vec.data() ptr ref, so this is a reproduction that needs to be kept updated

        Set set;
    };
}  // namespace cav

#endif

/* ################################################################# */
/* #### Original Header: include/UniqueCol.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_ROUTE_HPP_
#define SPH_INCLUDE_ROUTE_HPP_

#include <algorithm>

/* #include "IndexList.hpp" */
/* #include "cft.hpp" */

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

/* ################################################################# */
/* #### Original Header: include/UniqueColSet.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_ROUTEPOOL_HPP_
#define SPH_INCLUDE_ROUTEPOOL_HPP_

#include <cassert>

/* #include "UniqueCol.hpp" */
/* #include "VectorSet.hpp" */

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

/* ################################################################# */
/* #### Original Header: include/Instance.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_INSTANCE_HPP_
#define SPH_INCLUDE_INSTANCE_HPP_

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <cassert>
#include <numeric>

/* #include "CollectionOf.hpp" */
/* #include "MStar.hpp" */
/* #include "Stopwatch.hpp" */
/* #include "UniqueColSet.hpp" */
/* #include "cft.hpp" */

#define SUBINST_MIN_COV 4U
#define SUBINST_MIN_SOLCOST_COV 4U
#define SUBINST_HARD_CAP 15'000U

/**
 * @brief Represente a complete instance of a Set Partitioning problem.
 *
 */
class Instance {
public:
    explicit Instance(const idx_t nrows_) : nrows(nrows_), active_rows(nrows, true), nactive_rows(nrows), fixed_cost(0.0) { }

    [[nodiscard]] inline auto get_ncols() const { return cols.size(); }
    [[nodiscard]] inline idx_t get_nrows() const { return nrows; }
    [[nodiscard]] inline idx_t get_active_rows_size() const { return nactive_rows; }

    [[nodiscard]] inline auto &get_active_cols() { return active_cols; }
    [[nodiscard]] inline auto &get_fixed_cols() { return fixed_cols; }
    [[nodiscard]] inline auto &get_cols() { return cols; }
    [[nodiscard]] inline Column &get_col(idx_t idx) { return cols[idx]; }
    [[nodiscard]] inline const Column &get_col(idx_t idx) const { return cols[idx]; }
    [[nodiscard]] inline real_t get_fixed_cost() { return fixed_cost; }

    inline void set_timelimit(double seconds) { timelimit = Timer(seconds); }
    [[nodiscard]] inline Timer &get_timelimit() { return timelimit; }

    [[nodiscard]] inline bool is_row_active(idx_t gi) {
        assert(gi < nrows);
        return active_rows[gi];
    }

    void inline fix_columns(const std::vector<idx_t> &idxs) {
        for (idx_t j : idxs) {
            for (idx_t i : cols[j]) { active_rows[i] = false; }
        }

        _fix_columns(idxs);
    }

    void reset_fixing() {
        assert(std::is_sorted(active_cols.begin(), active_cols.end()));
        assert(std::is_sorted(fixed_cols.begin(), fixed_cols.end()));

        active_rows.assign(nrows, true);

        // merge active and fixed columns
        active_cols.resize(cols.size());
        std::iota(active_cols.begin(), active_cols.end(), 0);

        fixed_cols.clear();
        fixed_cost = 0.0;
    }

    void fix_columns(const std::vector<idx_t> &idxs, const MStar &M_star) {
        assert(fixed_cols.empty());
        assert(active_rows.size() == nrows);
        assert(M_star.size() == nrows);

        for (idx_t i = 0; i < nrows; ++i) { active_rows[i] = !M_star[i]; }
        nactive_rows = M_star.get_uncovered();

        _fix_columns(idxs);
    }

    template <typename... _Args>
    idx_t add_column(_Args &&...args) {
        if (cols.add_column(std::forward<_Args>(args)...)) {
            active_cols.emplace_back(cols.size() - 1);
            return cols.size() - 1;
        }
        return REMOVED_INDEX;
    }

    template <typename ColContainer>
    std::vector<idx_t> add_columns(const ColContainer &new_cols) {

        std::vector<idx_t> inserted_cols_idxs;
        inserted_cols_idxs.reserve(new_cols.size());

        for (auto &new_col : new_cols) {
            idx_t inserted_idx = add_column(new_col);
            if (inserted_idx != REMOVED_INDEX) { inserted_cols_idxs.emplace_back(inserted_idx); }
        }

        return inserted_cols_idxs;
    }

    template <typename ColContainer>
    std::vector<idx_t> add_columns(ColContainer &&new_cols) {

        std::vector<idx_t> inserted_cols_idxs;
        inserted_cols_idxs.reserve(new_cols.size());

        for (auto &new_col : new_cols) {
            idx_t inserted_idx = add_column(std::move(new_col));
            if (inserted_idx != REMOVED_INDEX) { inserted_cols_idxs.emplace_back(inserted_idx); }
        }

        return inserted_cols_idxs;
    }

    template <typename CostVec, typename SolCostVec, typename MatBegVec, typename MatValVec>
    std::vector<idx_t> add_columns(const CostVec &costs, const SolCostVec &sol_costs, const MatBegVec &matbeg, const MatValVec &matval) {
        assert(costs.size() == sol_costs.size() && costs.size() <= matbeg.size());

        idx_t ncols = costs.size();

        std::vector<idx_t> inserted_cols_idxs;
        inserted_cols_idxs.reserve(costs.size());

        for (idx_t j = 0; j < ncols - 1; ++j) {
            idx_t inserted_idx = add_column(&matval[matbeg[j]], &matval[matbeg[j + 1]], costs[j], sol_costs[j]);
            if (inserted_idx != REMOVED_INDEX) { inserted_cols_idxs.emplace_back(inserted_idx); }
        }

        // Matbeg might prbably contains only the starts
        idx_t inserted_idx = add_column(&matval[matbeg[ncols - 1]], &matval[matval.size() - 1], costs[ncols - 1], sol_costs[ncols - 1]);
        if (inserted_idx != REMOVED_INDEX) { inserted_cols_idxs.emplace_back(inserted_idx); }

        assert(cols.size() == costs.size());
        return inserted_cols_idxs;
    }

    void fill_with_best_columns(std::vector<idx_t> &global_idxs) {
        _init_priced_cols(priced_cols);
        covering_times.reset_uncovered(nrows);
        std::sort(priced_cols.begin(), priced_cols.end(), [](const Priced_Col &c1, const Priced_Col &c2) { return c1.c_u < c2.c_u; });

        _select_C2_cols(priced_cols, covering_times, global_idxs);
        _select_C3_cols(priced_cols, global_idxs);
    }

    real_t fill_with_best_columns(std::vector<idx_t> &global_idxs, const std::vector<real_t> &u_k) {

        real_t global_LB = _price_active_cols(u_k, priced_cols);

        covering_times.reset_uncovered(nrows);  // reset convered_rows to consider only reduced costs covering for C2
        std::sort(priced_cols.begin(), priced_cols.end(), [](const Priced_Col &c1, const Priced_Col &c2) { return c1.c_u < c2.c_u; });

        _select_C1_cols(priced_cols, covering_times, global_idxs);
        _select_C2_cols(priced_cols, covering_times, global_idxs);
        _select_C3_cols(priced_cols, global_idxs);

        return global_LB;
    }

    /**
     * @brief Prune columns maintaining only the best ones.
     *
     * @tparam Hard_cap
     * @param u_k
     * @return std::vector<idx_t> map from old indexes to new ones to translate pre-existing solutions.
     *          REMOVED_INDEX if the column has been removed.
     */
    template <unsigned long Hard_cap>
    std::vector<idx_t> prune_instance(const std::vector<real_t> &u_k) {
        if (cols.size() > 3 * Hard_cap) {
            std::vector<idx_t> idxs_to_keep;
            idxs_to_keep.reserve(Hard_cap);

            _price_active_cols(u_k, priced_cols);
            covering_times.reset_uncovered(nrows);  // reset convered_rows to consider only reduced costs covering for C2

            _select_C1_cols<MAX_INDEX, Hard_cap>(priced_cols, covering_times, idxs_to_keep);
            _select_C2_cols<MAX_INDEX, Hard_cap>(priced_cols, covering_times, idxs_to_keep);
            _select_C3_cols<MAX_INDEX, Hard_cap>(priced_cols, idxs_to_keep);

            UniqueColSet new_cols;
            new_cols.reserve(idxs_to_keep.size());
            std::vector<idx_t> old_to_new_idx_map(cols.size(), REMOVED_INDEX);

            for (idx_t gj : idxs_to_keep) {
                old_to_new_idx_map[gj] = new_cols.size();
                new_cols.add_column(cols[gj]);
            }
            std::swap(cols, new_cols);

            return old_to_new_idx_map;
        }

        return std::vector<idx_t>();
    }


private:
    struct Priced_Col {
        idx_t j;
        real_t c_u;
        real_t sol_cost;
    };

    class Priced_Columns : public std::vector<Priced_Col> {
    public:
        Priced_Columns() { }

        void reset(idx_t ncols) {
            assert(ncols > 0);
            resize(ncols);
        }

        inline void select(idx_t n) {
            (*this)[n].j = REMOVED_INDEX;
            (*this)[n].c_u = (*this)[n].sol_cost = REAL_MAX;
        }
        inline bool is_selected(idx_t n) const { return (*this)[n].j == REMOVED_INDEX; }
    };

    void _init_priced_cols(Priced_Columns &_priced_cols) {
        _priced_cols.reset(active_cols.size());

        idx_t p_idx = 0;
        for (idx_t gj : active_cols) {
            assert(gj < cols.size());

            Column &col = cols[gj];
            for (idx_t gi : col) {
                if (active_rows[gi]) {
                    _priced_cols[p_idx++] = {gj, col.get_cost(), col.get_solcost()};
                    break;
                }
            }
        }
        _priced_cols.resize(p_idx);
    }

    real_t _price_active_cols(const std::vector<real_t> &u_k, Priced_Columns &_priced_cols) {

        _priced_cols.reset(active_cols.size());

        // price all active columns and add their contribution to the LB
        real_t global_LB = std::reduce(u_k.begin(), u_k.end(), static_cast<real_t>(0.0));

        idx_t p_idx = 0;
        for (idx_t gj : active_cols) {
            assert(gj < cols.size());

            Column &col = cols[gj];
            real_t c_u = col.get_cost();

            bool is_empty = true;
            for (idx_t gi : col) {
                if (active_rows[gi]) {
                    is_empty = false;
                    c_u -= u_k[gi];  // NOTE: multipliers need to be adapted to global multipliers!!!!!
                }
            }

            if (!is_empty) {  // check for empty columns
                if (c_u < 0.0) { global_LB += c_u; }

                _priced_cols[p_idx++] = {gj, c_u, col.get_solcost()};
            }
        }

        _priced_cols.resize(p_idx);

        return global_LB;
    }

    template <unsigned long Min_cov = SUBINST_MIN_COV, unsigned long Hard_cap = SUBINST_HARD_CAP>
    NO_INLINE void _select_C1_cols(Priced_Columns &_priced_cols, MStar &_covering_times, std::vector<idx_t> &global_col_idxs) {

        idx_t fivem = std::min<idx_t>(Hard_cap, std::min<idx_t>(Min_cov * nactive_rows, _priced_cols.size()));
        global_col_idxs.reserve(fivem);

        std::sort(_priced_cols.begin(), _priced_cols.end(), [](const Priced_Col &c1, const Priced_Col &c2) { return c1.c_u < c2.c_u; });

        for (idx_t n = 0; n < fivem; n++) {
            assert(n < _priced_cols.size());

            if (_priced_cols.is_selected(n) || _priced_cols[n].c_u >= 0.1) { continue; }

            idx_t gj = _priced_cols[n].j;
            assert(gj < cols.size());
            assert(std::count_if(cols[gj].begin(), cols[gj].end(), [&](idx_t i) { return active_rows[i]; }) > 0);

            global_col_idxs.emplace_back(gj);
            _covering_times.cover_rows(cols[gj]);
            _priced_cols.select(n);
        }

        IF_DEBUG {
            [[maybe_unused]] auto old_end = global_col_idxs.end();
            assert(std::unique(global_col_idxs.begin(), global_col_idxs.end()) == old_end);
        }
    }

    template <unsigned long Min_cov = SUBINST_MIN_COV, unsigned long Hard_cap = SUBINST_HARD_CAP>
    NO_INLINE void _select_C2_cols(Priced_Columns &_priced_cols, MStar &_covering_times, std::vector<idx_t> &global_col_idxs) {

        assert(std::is_sorted(_priced_cols.begin() + global_col_idxs.size(), _priced_cols.end(),
                              [](const Priced_Col &c1, const Priced_Col &c2) { return c1.c_u < c2.c_u; }));

        if (nactive_rows == 0) { }

        idx_t min_cov = std::min<idx_t>(Min_cov, Hard_cap / nactive_rows);
        idx_t fivem = std::min<idx_t>(min_cov * nactive_rows, _priced_cols.size());
        global_col_idxs.reserve(fivem);

        // check for still-uncovered rows
        idx_t rows_to_cover = 0;
        for (idx_t gi = 0; gi < nrows; ++gi) {
            if (active_rows[gi]) {
                _covering_times[gi] = min_cov - std::min<idx_t>(min_cov, _covering_times[gi]);
                rows_to_cover += static_cast<idx_t>(_covering_times[gi] > 0);
            } else {
                _covering_times[gi] = 0;
            }
        }

        for (idx_t n = global_col_idxs.size(); n < _priced_cols.size(); ++n) {
            assert(!_priced_cols.is_selected(n));

            Column &col = cols[_priced_cols[n].j];
            for (idx_t gi : col) {
                if (_covering_times[gi] == 0) { continue; }

                --_covering_times[gi];

                if (!_priced_cols.is_selected(n)) {
                    idx_t gj = _priced_cols[n].j;
                    assert(gj < cols.size());
                    assert(std::count_if(cols[gj].begin(), cols[gj].end(), [&](idx_t i) { return active_rows[i]; }) > 0);

                    global_col_idxs.emplace_back(gj);
                    _priced_cols.select(n);
                }

                rows_to_cover -= static_cast<idx_t>(_covering_times[gi] == 0);
                if (rows_to_cover == 0) {
                    assert(std::count(_covering_times.begin(), _covering_times.end(), 0) == _covering_times.size());
                    break;
                }
            }
        }

        IF_DEBUG {
            [[maybe_unused]] auto old_end = global_col_idxs.end();
            assert(std::unique(global_col_idxs.begin(), global_col_idxs.end()) == old_end);
        }
    }

    template <unsigned long Min_cov = SUBINST_MIN_SOLCOST_COV, unsigned long Hard_cap = SUBINST_HARD_CAP>
    NO_INLINE void _select_C3_cols(Priced_Columns &_priced_cols, std::vector<idx_t> &global_col_idxs) {
        idx_t fivem = std::min<idx_t>(Hard_cap, std::min<idx_t>(Min_cov * nactive_rows, _priced_cols.size()));
        global_col_idxs.reserve(fivem);

        std::nth_element(_priced_cols.begin(), _priced_cols.begin() + fivem, _priced_cols.end(),
                         [](const Priced_Col &c1, const Priced_Col &c2) { return c1.sol_cost < c2.sol_cost; });

        if (_priced_cols[0].sol_cost == REAL_MAX) { return; }

        for (idx_t n = 0; n < fivem; ++n) {
            assert(n < _priced_cols.size());

            if (_priced_cols.is_selected(n) || _priced_cols[n].sol_cost == REAL_MAX) { continue; }

            idx_t gj = _priced_cols[n].j;
            assert(gj < cols.size());
            assert(std::count_if(cols[gj].begin(), cols[gj].end(), [&](idx_t i) { return active_rows[i]; }) > 0);

            global_col_idxs.emplace_back(gj);
            _priced_cols.select(n);
        }

        IF_DEBUG {
            [[maybe_unused]] auto old_end = global_col_idxs.end();
            assert(std::unique(global_col_idxs.begin(), global_col_idxs.end()) == old_end);
        }
    }

    void _fix_columns(const std::vector<idx_t> &idxs) {
        idx_t iok = 0;
        for (idx_t j = 0; j < cols.size(); ++j) {
            auto &col = cols[j];
            for (idx_t i : col) {
                if (active_rows[i]) {
                    active_cols[iok++] = j;
                    break;
                }
            }
        }

        active_cols.resize(iok);
        fixed_cols = idxs;

        fixed_cost = 0.0;
        for (idx_t j : fixed_cols) { fixed_cost += cols[j].get_cost(); }
    }


private:
    const idx_t nrows;
    UniqueColSet cols;
    std::vector<idx_t> active_cols;
    std::vector<idx_t> fixed_cols;
    std::vector<bool> active_rows;
    idx_t nactive_rows;
    real_t fixed_cost;
    Priced_Columns priced_cols;
    MStar covering_times;

    Timer timelimit;
};

#endif

/* ################################################################# */
/* #### Original Header: include/SubInstance.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_SUBINSTANCE_HPP_
#define SPH_INCLUDE_SUBINSTANCE_HPP_


#include <fmt/core.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <cassert>
#include <numeric>

/* #include "CollectionOf.hpp" */
/* #include "IndexList.hpp" */
/* #include "Instance.hpp" */
/* #include "MStar.hpp" */
/* #include "Stopwatch.hpp" */
/* #include "cft.hpp" */

class SubInstance {

public:
    explicit SubInstance(Instance &inst_) : inst(inst_), fixed_cost(inst_.get_fixed_cost()) { }

    [[nodiscard]] inline auto get_ncols() const { return cols.size(); }
    [[nodiscard]] inline auto get_nrows() const { return rows.size(); }

    [[nodiscard]] inline auto get_global_col_idx(idx_t local_j) const { return local_to_global_col_idxs[local_j]; }
    [[nodiscard]] inline auto get_global_row_idx(idx_t local_i) const { return local_to_global_row_idxs[local_i]; }
    [[nodiscard]] inline auto get_local_row_idx(idx_t global_i) const { return global_to_local_row_idxs[global_i]; }

    [[nodiscard]] inline auto &get_cols() { return cols; }
    [[nodiscard]] inline auto &get_rows() { return rows; }

    [[nodiscard]] inline const auto &get_cols() const { return cols; }
    [[nodiscard]] inline const auto &get_rows() const { return rows; }

    [[nodiscard]] inline auto &get_col(idx_t idx) { return cols[idx]; }
    [[nodiscard]] inline auto &get_row(idx_t idx) { return rows[idx]; }

    [[nodiscard]] inline const auto &get_col(idx_t idx) const { return cols[idx]; }
    [[nodiscard]] inline const auto &get_row(idx_t idx) const { return rows[idx]; }

    [[nodiscard]] inline auto &get_fixed_cols() { return fixed_cols_global_idxs; }
    [[nodiscard]] inline auto get_fixed_cost() const { return fixed_cost; }
    [[nodiscard]] inline auto &get_instance() { return inst; }

    [[nodiscard]] Timer &get_timelimit() { return inst.get_timelimit(); }

    [[nodiscard]] bool is_corrupted() const {
        idx_t j_counter = 0;
        for (auto &col : cols) {
            if (std::addressof(col) != std::addressof(cols[j_counter])) {
                fmt::print("Subinstance cols iterator corrupted at {}: {} != {}\n", j_counter, (void *)std::addressof(col),
                           (void *)std::addressof(cols[j_counter]));
                return true;
            }
            ++j_counter;
        }

        for (idx_t j = 0; j < cols.size(); ++j) {
            if (cols[j].empty()) {
                IF_DEBUG { fmt::print("Col {} is empty.\n ", j); }
                return true;
            }
            if (j > get_ncols()) {
                IF_DEBUG { fmt::print("Col {} does not exist. \n Col: ", j, fmt::join(cols[j], ", ")); }
                return true;
            }

            for (idx_t i : cols[j]) {
                if (std::find(rows[i].begin(), rows[i].end(), j) == rows[i].end()) {
                    IF_DEBUG { fmt::print("Col {} not found in row {}. \n Row: ", j, i, fmt::join(rows[i], ", ")); }
                    return true;
                }
            }
        }

        for (idx_t i = 0; i < rows.size(); ++i) {
            if (rows[i].empty()) {
                IF_DEBUG { fmt::print("Row {} is empty.\n ", i); }
                return true;
            }
            if (i > get_nrows()) {
                IF_DEBUG { fmt::print("Row {} does not exist. \n Row: ", i, fmt::join(rows[i], ", ")); }
                return true;
            }

            for (idx_t j : rows[i]) {
                if (std::find(cols[j].begin(), cols[j].end(), i) == cols[j].end()) {
                    IF_DEBUG { fmt::print("Row {} not found in col {}. \n Col: ", i, j, fmt::join(cols[j], ", ")); }
                    return true;
                }
            }
        }

        return false;
    }

    [[nodiscard]] real_t get_global_LB(const std::vector<real_t> &u_k) {
        real_t global_LB = std::reduce(u_k.begin(), u_k.end(), static_cast<real_t>(0.0));

        // price all active columns and add their contribution to the LB
        for (idx_t gj : inst.get_active_cols()) {
            const auto &col = inst.get_col(gj);
            real_t c_u = col.get_cost();

            for (idx_t gi : col) {
                if (_is_global_row_active(gi)) {
                    idx_t li = global_to_local_row_idxs[gi];  // retrieve the mapped row index
                    c_u -= u_k[li];
                }
            }

            if (c_u < 0.0) { global_LB += c_u; }
        }

        return global_LB;
    }
    [[nodiscard]] idx_t find_local_col_idx(idx_t gj) {

        for (idx_t gi : inst.get_col(gj)) {
            if (_is_global_row_active(gi)) {
                idx_t active_li = global_to_local_row_idxs[gi];
                for (idx_t lj : rows[active_li]) {
                    if (local_to_global_col_idxs[lj] == gj) { return lj; }
                }
                break;
            }
        }
        return REMOVED_INDEX;
    }
    [[nodiscard]] std::vector<idx_t> get_localized_solution(const std::vector<idx_t> &glob_sol) {
        assert(glob_sol.size() >= fixed_cols_global_idxs.size());

        std::vector<idx_t> local_sol;
        local_sol.reserve(glob_sol.size() - fixed_cols_global_idxs.size() - inst.get_fixed_cols().size());
        for (idx_t gj : glob_sol) {
            idx_t lj = find_local_col_idx(gj);
            if (lj != REMOVED_INDEX) { local_sol.emplace_back(lj); }
        }

        IF_DEBUG {
            auto &sifc = fixed_cols_global_idxs;
            auto &ifc = inst.get_fixed_cols();
            for ([[maybe_unused]] idx_t lj : local_sol) { assert(std::find(glob_sol.begin(), glob_sol.end(), local_to_global_col_idxs[lj]) != glob_sol.end()); }
            for ([[maybe_unused]] idx_t gj : sifc) { assert(std::find(glob_sol.begin(), glob_sol.end(), gj) != glob_sol.end()); }
            for ([[maybe_unused]] idx_t gj : ifc) { assert(std::find(glob_sol.begin(), glob_sol.end(), gj) != glob_sol.end()); }
            for ([[maybe_unused]] idx_t gj : glob_sol) {

                [[maybe_unused]] bool check1 = [&]() {
                    for (idx_t lj : local_sol)
                        if (local_to_global_col_idxs[lj] == gj) return true;
                    return false;
                }();
                [[maybe_unused]] bool check2 = std::find(sifc.begin(), sifc.end(), gj) != sifc.end();
                [[maybe_unused]] bool check3 = std::find(ifc.begin(), ifc.end(), gj) != ifc.end();

                assert(check1 || check2 || check3);
            }
            assert(sifc.size() + ifc.size() + local_sol.size() == glob_sol.size());
        }

        return local_sol;
    }

    template <typename ColContainer>
    std::vector<idx_t> add_columns(const ColContainer &new_cols) {

        std::vector<idx_t> inserted_cols_idxs = inst.add_columns(new_cols);
        assert(!inserted_cols_idxs.empty());

        idx_t new_ncols = cols.size() + inserted_cols_idxs.size();
        idx_t old_ncols = cols.size();

        cols.reserve(new_ncols);
        local_to_global_col_idxs.resize(new_ncols);

        idx_t lj = old_ncols;
        for (idx_t &gj : inserted_cols_idxs) {

            local_to_global_col_idxs[lj] = gj;

            Column &gcol = inst.get_col(gj);
            cols.new_col_create(gcol.get_cost());
            for (idx_t gi : gcol) {
                assert(_is_global_row_active(gi));
                idx_t li = global_to_local_row_idxs[gi];
                cols.new_col_push_back(li);
                rows[li].emplace_back(lj);
            }

            gj = lj;  // convert global to local indices
            ++lj;
        }
        assert(lj == new_ncols);
        assert(!is_corrupted());

        return inserted_cols_idxs;
    }

    inline void update_sol_cost(const std::vector<idx_t> &local_sol) {
        real_t sol_cost = std::accumulate(local_sol.begin(), local_sol.end(), 0.0, [&](real_t sum, idx_t lj) { return sum + cols[lj].get_cost(); });
        sol_cost += get_fixed_cost();
        update_sol_costs(local_sol, sol_cost);
    }

    void update_sol_costs(const std::vector<idx_t> &local_sol, real_t sol_cost) {

        for (idx_t lj : local_sol) {
            Column &gcol = inst.get_col(local_to_global_col_idxs[lj]);
            if (gcol.get_solcost() > sol_cost) { gcol.set_solcost(sol_cost); }
        }

        for (idx_t gj : fixed_cols_global_idxs) {
            if (inst.get_col(gj).get_solcost() > sol_cost) { inst.get_col(gj).set_solcost(sol_cost); }
        }

        for (idx_t gj : inst.get_fixed_cols()) {
            if (inst.get_col(gj).get_solcost() > sol_cost) { inst.get_col(gj).set_solcost(sol_cost); }
        }
    }

    void reset() {
        local_to_global_row_idxs.resize(inst.get_nrows());
        global_to_local_row_idxs.resize(inst.get_nrows());

        idx_t li = 0;
        for (idx_t gi = 0; gi < inst.get_nrows(); ++gi) {
            if (inst.is_row_active(gi)) {
                local_to_global_row_idxs[li] = gi;
                global_to_local_row_idxs[gi] = li;
                ++li;
            } else {
                global_to_local_row_idxs[gi] = REMOVED_INDEX;
            }
        }
        local_to_global_row_idxs.resize(li);

        local_to_global_col_idxs.clear();
        fixed_cols_global_idxs.clear();
        fixed_cost = inst.get_fixed_cost();

        inst.fill_with_best_columns(local_to_global_col_idxs);
        replace_columns(local_to_global_col_idxs);

        IF_VERBOSE { fmt::print("Sub-instance size = {}x{}.\n", rows.size(), cols.size()); }

        assert(!is_corrupted());
    }

    NO_INLINE real_t price(const std::vector<real_t> &u_k) {

        global_u_k.assign(inst.get_nrows(), 0);
        for (idx_t i = 0; i < u_k.size(); ++i) { global_u_k[local_to_global_row_idxs[i]] = u_k[i]; }

        local_to_global_col_idxs.clear();
        real_t global_LB = inst.fill_with_best_columns(local_to_global_col_idxs, global_u_k);

        replace_columns(local_to_global_col_idxs);

        return global_LB;
    }

    idx_t fix_columns(std::vector<idx_t> &local_idxs_to_fix, std::vector<real_t> &u_star) {

        if (local_idxs_to_fix.empty()) { return rows.size(); }

        // mark rows to remove
        for (idx_t lj : local_idxs_to_fix) {
            idx_t gj = local_to_global_col_idxs[lj];
            local_to_global_col_idxs[lj] = REMOVED_INDEX;

            assert(gj != REMOVED_INDEX);
            assert(std::find(fixed_cols_global_idxs.begin(), fixed_cols_global_idxs.end(), gj) == fixed_cols_global_idxs.end());

            fixed_cols_global_idxs.emplace_back(gj);
            const auto &col = cols[lj];
            for (idx_t li : col) {
                if (local_to_global_row_idxs[li] == REMOVED_INDEX) { continue; }

                idx_t gi = local_to_global_row_idxs[li];
                local_to_global_row_idxs[li] = REMOVED_INDEX;
                global_to_local_row_idxs[gi] = REMOVED_INDEX;
            }
        }

        fixed_cost = inst.get_fixed_cost();
        for (idx_t gj : fixed_cols_global_idxs) { fixed_cost += inst.get_col(gj).get_cost(); }

        // compact rows
        idx_t li = 0;
        while (local_to_global_row_idxs[li] != REMOVED_INDEX) { ++li; }

        idx_t rows_left = li;
        for (; li < rows.size(); ++li) {
            if (local_to_global_row_idxs[li] == REMOVED_INDEX) { continue; }

            idx_t gi = local_to_global_row_idxs[li];

            assert(_is_global_row_active(gi));
            assert(!rows[rows_left].empty());

            global_to_local_row_idxs[gi] = rows_left;
            local_to_global_row_idxs[rows_left] = gi;
            u_star[rows_left] = u_star[li];
            ++rows_left;
        }
        local_to_global_row_idxs.resize(rows_left);
        u_star.resize(rows_left);

        IF_DEBUG {
            for ([[maybe_unused]] idx_t gi : local_to_global_row_idxs) { assert(_is_global_row_active(gi)); }
        }

        if (rows_left == 0) {
            cols.clear();
            rows.clear();
            return 0;
        }

        idx_t lj = 0;
        for (idx_t gj : local_to_global_col_idxs) {
            if (gj == REMOVED_INDEX) { continue; }

            for (auto gi : inst.get_col(gj)) {
                if (_is_global_row_active(gi)) {
                    local_to_global_col_idxs[lj++] = gj;
                    break;
                }
            }
        }
        local_to_global_col_idxs.resize(lj);

        replace_columns(local_to_global_col_idxs);
        return rows_left;
    }

    void replace_columns(const std::vector<idx_t> &glob_cols_idxs) {
        assert(!glob_cols_idxs.empty());

        rows.resize(local_to_global_row_idxs.size());
        for (auto &row : rows) { row.clear(); }

        idx_t ncols = glob_cols_idxs.size();

        cols.clear();
        cols.reserve(ncols);

        for (idx_t lj = 0; lj < glob_cols_idxs.size(); ++lj) {
            idx_t gj = glob_cols_idxs[lj];

            assert(std::find(fixed_cols_global_idxs.begin(), fixed_cols_global_idxs.end(), gj) == fixed_cols_global_idxs.end());
            assert(!inst.get_col(gj).empty());

            const auto &gcol = inst.get_col(gj);

            cols.new_col_create(gcol.get_cost());
            for (idx_t gi : gcol) {
                if (_is_global_row_active(gi)) {
                    idx_t li = global_to_local_row_idxs[gi];
                    cols.new_col_push_back(li);
                    rows[li].emplace_back(lj);
                }
            }
            assert(!cols[lj].empty());
        }
        assert(!is_corrupted());
    }


private:
    [[nodiscard]] inline bool _is_global_row_active(const idx_t gbl_idx) {
        assert(gbl_idx < global_to_local_row_idxs.size());
        return global_to_local_row_idxs[gbl_idx] != REMOVED_INDEX;
    }

private:
    Instance &inst;
    SubInstCols cols;
    std::vector<Row> rows;
    std::vector<idx_t> local_to_global_col_idxs;  // map local to original indexes
    std::vector<idx_t> fixed_cols_global_idxs;    // original indexes of locally fixed cols

    std::vector<idx_t> global_to_local_row_idxs;
    std::vector<idx_t> local_to_global_row_idxs;

    std::vector<real_t> global_u_k;
    real_t fixed_cost;
};

#endif  // SPH_INCLUDE_SUBINSTANCE_HPP_


/* ################################################################# */
/* #### Original Header: include/Solution.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_SOLUTION_HPP_
#define SPH_INCLUDE_SOLUTION_HPP_
#include <vector>

/* #include "SubInstance.hpp" */
/* #include "cft.hpp" */

/**
 * Solution classes are simply std::vector of indexes.
 * In order to distinguish between solutions belonging different parts of the algorithm, solution of different "type" are actually
 * different classes (in this way we can always know what a solution represent exactly).
 * However, in the end they are simply a vector of indexes plus some utility functions attached to them.
 */
class LocalSolution;
class LocalSolution;
class GlobalSolution;

class BaseSolution : public std::vector<idx_t> {
public:
    template <typename... Args>
    explicit BaseSolution(Args&&... _args) : std::vector<idx_t>(std::forward<Args>(_args)...) { }

    inline void mark_removal(idx_t idx) { (*this)[idx] = REMOVED_INDEX; }

    inline void apply_removal() {
        idx_t j2 = 0;
        for (idx_t j1 = 0; j1 < this->size(); ++j1) {
            if ((*this)[j1] != REMOVED_INDEX) { (*this)[j2++] = (*this)[j1]; }
        }
        this->resize(j2);
    }

    inline void remove(std::vector<idx_t>& col_idxs) {
        // n: sol size, n': cols to remove, assuming n >> n'

        std::sort(col_idxs.begin(), col_idxs.end());  // O(n' log(n'))
        idx_t removed_counter = col_idxs.size();

        idx_t j1 = 0, j2 = 0;
        while (j2 < this->size()) {  // O(n log(n'))
            if (std::binary_search(col_idxs.begin(), col_idxs.end(), (*this)[j2])) {
                (*this)[j2] = this->back();
                this->pop_back();
                if (--removed_counter == 0) { break; }

            } else {
                (*this)[j1++] = (*this)[j2++];
            }
        }

        while (j2 < this->size()) { (*this)[j1++] = (*this)[j2++]; }  // O(n)
        this->resize(j1);
    }

protected:
    template <typename Inst>
    [[nodiscard]] inline real_t compute_cost(Inst& inst) const {
        real_t cost = 0.0;
        for (auto j : *this) { cost += inst.get_col(j).get_cost(); }
        return cost;
    }
};

/**
 * @brief Solution used inside the 3-phase cycle. Use local indexes mappings (defined in the current subinstance) and does not
 * contain fixed columns.
 *
 */
class LocalSolution : public BaseSolution {
public:
    template <typename... Args>
    explicit LocalSolution(Args&&... _args) : BaseSolution(std::forward<Args>(_args)...) { }

    [[nodiscard]] inline real_t compute_cost(SubInstance& subinst) const { return BaseSolution::compute_cost(subinst); }
};


/**
 * @brief Complete solution coherent with the complete instance, can be used for I/O.
 *
 */
class GlobalSolution : public BaseSolution {
public:
    GlobalSolution() : cost(REAL_MAX) { }
    explicit GlobalSolution(idx_t size) : BaseSolution(size), cost(REAL_MAX) { }
    GlobalSolution(idx_t size, real_t val) : BaseSolution(size, val), cost(REAL_MAX) { }
    GlobalSolution(const GlobalSolution& other) : BaseSolution(other), cost(other.cost) { }
    GlobalSolution(GlobalSolution&& other)  noexcept : BaseSolution(std::move(other)), cost(other.cost) { }

    GlobalSolution& operator=(const GlobalSolution& other) = default;
    GlobalSolution& operator=(GlobalSolution&& other) = default;

    GlobalSolution(SubInstance& subinst, const LocalSolution& sol) {
        for (idx_t j : sol) { emplace_back(subinst.get_global_col_idx(j)); }
        for (idx_t j : subinst.get_fixed_cols()) { emplace_back(j); }
        for (idx_t j : subinst.get_instance().get_fixed_cols()) { emplace_back(j); }

        cost = BaseSolution::compute_cost(subinst.get_instance());
    }

    [[nodiscard]] inline real_t get_cost() const { return cost; }

    void set_cost(real_t cost_) { cost = cost_; }

private:
    real_t cost;
};


#endif

/* ################################################################# */
/* #### Original Header: include/ExactSolver.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_EXACTSOLVER_HPP_
#define SPH_INCLUDE_EXACTSOLVER_HPP_

#include <fmt/core.h>
#include <fmt/ranges.h>
#include <ilcplex/cplex.h>

/* #include "Solution.hpp" */
/* #include "SubInstance.hpp" */
/* #include "cft.hpp" */

#define RESIZE_UP(vec, sz) ;    if (vec.size() < sz) { vec.resize(sz); }

#define ASSIGN_UP(vec, sz, val) ;    if (vec.size() < sz) { vec.assign(sz, val); }

#define SET_INT(P, VAL)                                                            ;    if (int res = 0; (res = CPXsetintparam(env, P, VAL))) {                        ;        fmt::print(stderr, "Error while setting " #P " parameter at  {} \n", VAL); ;        return res;                                                                ;    }

#define SET_DBL(P, VAL)                                                           ;    if (int res = 0; (res = CPXsetdblparam(env, P, VAL))) {                       ;        fmt::print(stderr, "Error while setting " #P " parameter at {} \n", VAL); ;        return res;                                                               ;    }

class ExactSolver {
public:
    ExactSolver() : env(CPXopenCPLEX(nullptr)) { }

    ~ExactSolver() { CPXcloseCPLEX(&env); }

    LocalSolution build_and_opt(SubInstance& subinst, LocalSolution& warmstart, Timer& time_limit) {

        lp = CPXcreateprob(env, nullptr, "exact");
        int res = 0;

        if ((res = build_model(subinst))) {
            fmt::print(stderr, "Error while building the model (errno: {})\n", res);
            return LocalSolution();
        }

        if ((res = set_warmstart(warmstart))) { fmt::print(stderr, "Error while setting warmstart (errno: {})\n", res); }

        if ((res = set_CPX_params(time_limit.seconds_until_end()))) {
            fmt::print(stderr, "Error while setting parameter (errno: {})\n", res);
            return LocalSolution();
        }

        IF_DEBUG {
            if ((res = CPXwriteprob(env, lp, "model.lp", nullptr))) { fmt::print(stderr, "Error while writing problem file(errno: {})\n", res); }
        }

        if ((res = CPXmipopt(env, lp))) {
            fmt::print("Cplex finished with error {}\n", res);
            return LocalSolution();
        }

        int stat = 0;
        double obj = 0.0;
        int ncols = subinst.get_ncols();
        RESIZE_UP(dbl_vals, ncols + 1U);
        LocalSolution sol;

        if (((res = CPXsolution(env, lp, &stat, &obj, dbl_vals.data(), nullptr, nullptr, nullptr) == 0))) {
            for (int i = 0; i < ncols; ++i) {
                if (dbl_vals[i] > 0.5) { sol.emplace_back(i); }
            }

            MStar coverage(subinst.get_nrows());
            coverage.reset_covered(subinst.get_cols(), sol, subinst.get_nrows());
            if (coverage.get_uncovered() > 0) {
                fmt::print(stderr, "Error, solution does not cover all the rows!\n");
                fmt::print(stderr, " Row coverage:\n{}\n", fmt::join(coverage, ", "));
                fmt::print(stderr, " Cols value: {}\n", fmt::join(dbl_vals, ", "));
                return LocalSolution();
            }
        }

        CPXfreeprob(env, &lp);

        return sol;
    }

private:
    int build_model(SubInstance& subinst) {
        int res;

        if ((res = CPXchgobjoffset(env, lp, subinst.get_fixed_cost()))) {
            fmt::print(stderr, "Error while setting obj func constant term! (errno: {})\n", res);
            return res;
        }

        if ((res = add_variables(subinst))) {
            fmt::print(stderr, "Error while creating columns! (errno: {})\n", res);
            return res;
        }

        if ((res = add_constraints(subinst))) {
            fmt::print(stderr, "Error creating rows! (errno: {})\n", res);
            return res;
        }

        return 0;
    }

    int add_variables(SubInstance& subinst) {
        SubInstCols& cols = subinst.get_cols();
        idx_t ncols = subinst.get_ncols();

        ASSIGN_UP(ctype, ncols, 'B');
        ASSIGN_UP(ones, ncols, 1.0);
        ASSIGN_UP(lb, ncols, 0.0);
        RESIZE_UP(dbl_vals, ncols);

        std::transform(cols.begin(), cols.end(), dbl_vals.begin(), [](auto& c) { return c.get_cost(); });

        return CPXnewcols(env, lp, ncols, dbl_vals.data(), lb.data(), ones.data(), ctype.data(), nullptr);
    }

    int add_constraints(SubInstance& subinst) {
        std::vector<Row>& rows = subinst.get_rows();
        idx_t nrows = subinst.get_nrows();

        RESIZE_UP(rmatbeg, nrows);
        rmatind.clear();
        rmatind.reserve(nrows * 5U);

        int nzcount = 0;
        for (idx_t i = 0; i < nrows; ++i) {
            rmatbeg[i] = nzcount;
            nzcount += rows[i].size();
            rmatind.insert(rmatind.end(), rows[i].begin(), rows[i].end());
        }

        ASSIGN_UP(ones, static_cast<size_t>(nzcount), 1.0);
        ASSIGN_UP(sense, static_cast<size_t>(nrows), 'E');

        return CPXaddrows(env, lp, 0, nrows, nzcount, ones.data(), sense.data(), rmatbeg.data(), rmatind.data(), ones.data(), nullptr, nullptr);
    }

    int set_warmstart(LocalSolution& warmstart) {
        idx_t wsize = warmstart.size();
        if (wsize > 0) {

            int zero_int = 0;
            int effort = CPX_MIPSTART_NOCHECK;
            ASSIGN_UP(ones, wsize, 1.0);
            RESIZE_UP(rmatind, wsize);
            for (idx_t n = 0; n < wsize; ++n) { rmatind[n] = warmstart[n]; }

            return CPXaddmipstarts(env, lp, 1, wsize, &zero_int, rmatind.data(), ones.data(), &effort, nullptr);
        }

        return 0;
    }

    int set_CPX_params(double tlim) {

        if (int res = 0; (res = CPXchgobjsen(env, lp, CPX_MIN))) {
            fmt::print(stderr, "Error while setting objective function ctype! (errno: {})\n", res);
            return res;
        }

        SET_INT(CPXPARAM_ScreenOutput, CPX_OFF);
        SET_INT(CPXPARAM_MIP_Display, 2);
        SET_INT(CPXPARAM_Threads, 1);
        SET_INT(CPXPARAM_Parallel, CPX_PARALLEL_OPPORTUNISTIC);
        SET_DBL(CPXPARAM_MIP_PolishAfter_Time, tlim * 0.5);
        SET_DBL(CPXPARAM_TimeLimit, tlim);
        SET_INT(CPXPARAM_Emphasis_MIP, CPX_MIPEMPHASIS_OPTIMALITY);

        IF_VERBOSE {
            IF_DEBUG { SET_INT(CPXPARAM_ScreenOutput, CPX_ON); }
        }

        return 0;
    }

    ////////// PRIVATE FIELDS //////////

    CPXENVptr env;
    CPXLPptr lp;

    // Avoid repeated memory allocations
    std::vector<double> lb;        // variables lower bounds (ncols, 0.0)
    std::vector<double> ones;      // constaints coefficients (nzcount, 1.0)
    std::vector<char> ctype;       // binary variables (ncols, 'B')
    std::vector<int> rmatind;      // constraints indices (nzcount)
    std::vector<int> rmatbeg;      // constraints begin in rmatind or rmatval (nrows)
    std::vector<char> sense;       // constraints sense (nrows, 'E')
    std::vector<double> dbl_vals;  // store variables final values (ncols)
};


#endif

/* ################################################################# */
/* #### Original Header: include/Multipliers.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_MULTIPLIERS_HPP_
#define SPH_INCLUDE_MULTIPLIERS_HPP_
#include <vector>

/* #include "Instance.hpp" */
/* #include "cft.hpp" */

class LocalMultipliers : public std::vector<real_t> {
public:
    template <typename... Args>
    explicit LocalMultipliers(Args&&... _args) : std::vector<real_t>(std::forward<Args>(_args)...) { }

    [[nodiscard]] inline real_t compute_lb(SubInstance& subinst) const {
        auto& cols = subinst.get_cols();
        real_t LB = std::reduce(begin(), end(), 0.0);
        for (auto& col : cols) {
            real_t c_u = col.compute_lagr_cost(*this);
            if (c_u < 0.0) { LB += c_u; }
        }

        return LB;
    }
};


class GlobalMultipliers : public std::vector<real_t> {
public:
    GlobalMultipliers() : lb(REAL_LOWEST) { }
    explicit GlobalMultipliers(idx_t size) : std::vector<real_t>(size), lb(REAL_LOWEST) { }
    GlobalMultipliers(idx_t size, real_t val) : std::vector<real_t>(size, val), lb(REAL_LOWEST) { }
    GlobalMultipliers(const GlobalMultipliers& other) : std::vector<real_t>(other), lb(other.lb) { }
    GlobalMultipliers(GlobalMultipliers&& other) noexcept : std::vector<real_t>(std::move(other)), lb(other.lb) { }

    GlobalMultipliers& operator=(const GlobalMultipliers& other) = default;
    GlobalMultipliers& operator=(GlobalMultipliers&& other) = default;

    GlobalMultipliers(SubInstance& subinst, LocalMultipliers mult) : std::vector<real_t>(subinst.get_instance().get_nrows(), 0.0), lb(REAL_LOWEST) {
        for (idx_t i = 0; i < mult.size(); ++i) { (*this)[subinst.get_global_row_idx(i)] = mult[i]; }
        update_lb(subinst.get_instance());
    }

    inline void update_lb(Instance& inst) {

        real_t lb1 = inst.get_fixed_cost();
        real_t lb2 = 0.0;
        for (idx_t gi = 0; gi < size(); ++gi) {
            if (inst.is_row_active(gi)) { lb2 += (*this)[gi]; }
        }

        real_t lb3 = 0.0;
        auto& cols = inst.get_cols();
        for (idx_t gj : inst.get_active_cols()) {

            real_t c_u = cols[gj].get_cost();
            for (idx_t gi : cols[gj]) {
                if (inst.is_row_active(gi)) { c_u -= (*this)[gi]; }
            }

            if (c_u < 0.0) { lb3 += c_u; }
        }

        lb = lb1 + lb2 + lb3;
        // fmt::print("fixed {}, u {}, c_u {} = lb {}\n", lb1, lb2, lb3, lb);
    }

    [[nodiscard]] inline real_t get_lb() const { return lb; }


private:
    real_t lb;
};


#endif

/* ################################################################# */
/* #### Original Header: include/LowerBound.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_LOWERBOUND_HPP_
#define SPH_INCLUDE_LOWERBOUND_HPP_

#include <algorithm>
#include <array>
#include <cassert>
#include <vector>

/* #include "MStar.hpp" */
/* #include "Multipliers.hpp" */
/* #include "Solution.hpp" */
/* #include "SubInstance.hpp" */
/* #include "cft.hpp" */

/**
 * @brief Return a LB cost computed by row.
 * For each row its "cost" is computed as the cost of the smallest column divided by the size of the column.
 *
 */
template <typename Iter>
real_t best_col_LB(SubInstance& subinst, MStar& covered_rows, const Iter& begin, const Iter& end) {

    auto& cols = subinst.get_cols();
    real_t LB = 0.0;

    std::vector<std::pair<idx_t, idx_t>> covered_rows_backup;
    for (auto it = begin; it != end; ++it) {
        idx_t covered_counter = 0;
        auto& col = cols[*it];
        for (auto i : col) {
            if (covered_rows[i] > 0) {
                covered_rows_backup.emplace_back(i, covered_rows[i]);
                covered_rows[i] = 0;
                ++covered_counter;
            }
        }

        LB += (col.get_cost() * covered_counter) / static_cast<real_t>(col.size());
    }

    for (auto& [i, bk] : covered_rows_backup) { covered_rows[i] = bk; }

    return LB;
}

template <typename Iter>
real_t best_col_LB(MStar& covered_rows, const Iter& begin, const Iter& end) {

    real_t LB = 0.0;

    std::vector<std::pair<idx_t, idx_t>> covered_rows_backup;
    for (auto it = begin; it != end; ++it) {
        idx_t covered_counter = 0;
        auto& col = *it;
        for (auto i : col) {
            if (covered_rows[i] > 0) {
                covered_rows_backup.emplace_back(i, covered_rows[i]);
                covered_rows[i] = 0;
                ++covered_counter;
            }
        }

        LB += (col.get_cost() * covered_counter) / static_cast<real_t>(col.size());
    }

    for (auto& [i, bk] : covered_rows_backup) { covered_rows[i] = bk; }

    return LB;
}

/**
 * @brief Return the LB computed using the Lagrangian multipliers.
 *
 */
template <typename Iter, int n = 0, int d = 1>
real_t lagr_mul_LB(SubInstance& subinst, const LocalMultipliers& u_k, const Iter& begin, const Iter& end) {
    constexpr real_t threshold = static_cast<real_t>(n) / static_cast<real_t>(d);

    real_t LB = std::reduce(u_k.begin(), u_k.end(), 0.0);

    for (auto it = begin; it != end; ++it) {
        auto& col = subinst.get_col(*it);
        real_t c_u = col.compute_lagr_cost(u_k);
        if (c_u < threshold) { LB += c_u; }
    }

    return LB;
}

template <int n = 0, int d = 1>  // seriously c++?
real_t lagr_mul_LB(SubInstance& subinst, const std::vector<real_t>& u_k) {
    constexpr real_t threshold = static_cast<real_t>(n) / static_cast<real_t>(d);

    auto& cols = subinst.get_cols();
    real_t LB = std::reduce(u_k.begin(), u_k.end(), 0.0);
    // fmt::print("BASE LB ORIG: {}\n", LB);
    for (auto& col : cols) {
        real_t c_u = col.compute_lagr_cost(u_k);
        if (c_u < threshold) { LB += c_u; }
    }

    return LB;
}

template <int n = 0, int d = 1>  // seriously c++?
real_t lagr_mul_LB(const SubInstance& subinst, const std::vector<real_t>& u_k) {
    constexpr real_t threshold = static_cast<real_t>(n) / static_cast<real_t>(d);

    const auto& cols = subinst.get_cols();
    real_t LB = std::reduce(u_k.begin(), u_k.end(), 0.0);
    // fmt::print("BASE LB ORIG: {}\n", LB);
    for (const auto& col : cols) {
        real_t c_u = col.compute_lagr_cost(u_k);
        if (c_u < threshold) { LB += c_u; }
    }

    return LB;
}

template <int n = 0, int d = 1>  // seriously c++?
class DeltaLowerBound {
    static constexpr real_t threshold = static_cast<real_t>(n) / static_cast<real_t>(d);

public:
    DeltaLowerBound() : LB(0.0) { }

    real_t compute(SubInstance& subinst, const std::vector<real_t>& u_k) { return LB = lagr_mul_LB<n, d>(subinst, u_k); }

    real_t update(SubInstance& subinst, const std::vector<std::pair<idx_t, real_t>>& delta_u_k) {

        auto& cols = subinst.get_cols();
        auto& rows = subinst.get_rows();
        // real_t base = 0.0;
        for (auto& [i, delta_ui] : delta_u_k) {
            LB += delta_ui;
            // base += delta_ui;
            for (auto j : rows[i]) {
                real_t old_c = cols[j].get_cu();
                real_t new_c = cols[j].update_cu(delta_ui);

                if (new_c < threshold) { LB += new_c; }
                if (old_c < threshold) { LB -= old_c; }
            }
        }

        // fmt::print("BASE LB DELT: {}\n", base);
        return LB;
    }

private:
    real_t LB;
};

template <int n = 1, int d = 1000>  // seriously c++?
class ReducedLagrMultLB {
    static constexpr real_t threshold = static_cast<real_t>(n) / static_cast<real_t>(d);

public:
    void compute_reduced_sol(SubInstance& subinst, LocalSolution& S, MStar& I_S) {
        const auto& cols = subinst.get_cols();
        idx_t ncols = cols.size();

        I_S.reset_uncovered(subinst.get_nrows());

        S.clear();
        for (idx_t j = 0; j < ncols; ++j) {
            if (cols[j].get_cu() <= threshold) {
                S.emplace_back(j);
                I_S.cover_rows(cols[j]);
            }
        }

        R.clear();
        idx_t Ssize = S.size();
        for (idx_t sj = 0; sj < Ssize; ++sj) {
            if (I_S.is_redundant(cols[S[sj]])) { R.emplace_back(sj); }
        }

        std::sort(R.begin(), R.end(), [&](const auto sj1, const auto sj2) { return cols[S[sj1]].get_cu() > cols[S[sj2]].get_cu(); });

        for (const auto j : R) {
            if (I_S.is_redundant(cols[S[j]])) {
                I_S.uncover_rows(cols[S[j]]);
                S.mark_removal(j);
            }
        }
        S.apply_removal();
    }

    void compute_sol(SubInstance& subinst, LocalSolution& S, MStar& I_S) {
        const auto& cols = subinst.get_cols();
        idx_t ncols = cols.size();

        S.clear();
        for (idx_t j = 0; j < ncols; ++j) {
            if (cols[j].get_cu() <= threshold) { S.emplace_back(j); }
        }

        I_S.reset_covered(cols, S, subinst.get_nrows());
    }

private:
    std::vector<idx_t> R;
};

#endif

/* ################################################################# */
/* #### Original Header: include/Counter.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_COUNTER_HPP_
#define SPH_INCLUDE_COUNTER_HPP_

/* #include "cft.hpp" */

class Counter {

public:

    explicit Counter(idx_t max) : m(max), i(0) { }

    inline void inc() {
        i++;
    }

    [[nodiscard]] inline auto reached() const {
        return i >= m;
    }

    inline void restart() {
        i = 0;
    }

    inline void set_max(idx_t new_max) {
        m = new_max;
    }

    [[nodiscard]] inline idx_t get_max() const {
        return m;
    }

    [[nodiscard]] inline auto get() const {
        return i;
    }

private:

    idx_t m; // max value
    idx_t i; // current value

};


#endif  // SPH_INCLUDE_COUNTER_HPP_


/* ################################################################# */
/* #### Original Header: include/SubGradientUtils.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_SUBGRADIENTUTILS_HPP_
#define SPH_INCLUDE_SUBGRADIENTUTILS_HPP_

/* #include "Counter.hpp" */
/* #include "SubInstance.hpp" */
/* #include "cft.hpp" */
class StepSizeFactor {
public:
    StepSizeFactor(real_t initial_value_ = 0.1, idx_t period = 20) : p(period), lambda(initial_value_),
                                                                     wrst_LB(REAL_MAX), best_LB(REAL_LOWEST) { }

    inline void update(real_t current_LB) {
        if (current_LB < wrst_LB) {
            wrst_LB = current_LB;
        }
        if (current_LB > best_LB) {
            best_LB = current_LB;
        }

        p.inc();
        if (p.reached()) {
            const auto perc_diff = (best_LB - wrst_LB) / best_LB;
            //assert(perc_diff > 0.0);

            if (perc_diff > 0.01) {
                lambda /= 2.0;
            } else if (perc_diff <= 0.001) {
                lambda *= 1.5;
            }

            p.restart();
            best_LB = REAL_LOWEST;
            wrst_LB = REAL_MAX;
        }
    }

    [[nodiscard]] inline real_t get() const { return lambda; }

    inline void reset() {
        p.restart();
        best_LB = REAL_LOWEST;
        wrst_LB = REAL_MAX;
    }

private:
    Counter p;
    real_t lambda;
    real_t wrst_LB;
    real_t best_LB;
};

class PricingPeriod : public Counter {
public:
    PricingPeriod(idx_t period = 10U, idx_t period_lim = 1000UL) : Counter(period), Tlim(period_lim) { }

    void restart() = delete;

    inline void reset(real_t global_LB, real_t local_LB, real_t UB) {
        Counter::restart();

        const auto delta = (local_LB - global_LB) / UB;

        if (delta <= 0.000001) {
            set_max(std::min<idx_t>(Tlim, get_max() * 10));
        } else if (delta <= 0.02) {
            set_max(std::min<idx_t>(Tlim, get_max() * 5));
        } else if (delta <= 0.2) {
            set_max(std::min<idx_t>(Tlim, get_max() * 2));
        } else {
            set_max(10);
        }
    }

private:
    idx_t Tlim;
};

class ExitCondition {
public:
    explicit ExitCondition(idx_t period_ = 300) : period(period_), counter(period_), LB_past(REAL_LOWEST) { }

    bool operator()(real_t LB) {
        if (--counter == 0) {
            const auto impr = LB - LB_past;
            const auto gap = 2.0 * (LB - LB_past) / LB;

            LB_past = LB;
            counter = period;

            return impr < 1 && gap < 0.001;
        }
        return false;
    }

private:
    const idx_t period;
    idx_t counter;
    real_t LB_past;
};

#endif

/* ################################################################# */
/* #### Original Header: include/SubGradient.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_SUBGRADIENT_HPP_
#define SPH_INCLUDE_SUBGRADIENT_HPP_

#include <cassert>
#include <random>

/* #include "LowerBound.hpp" */
/* #include "SubGradientUtils.hpp" */
/* #include "cft.hpp" */

#define REAL_TOLERANCE 1E-6

class SubGradient {

public:
    /**
     * Initialize a subgradient object by setting the initial multipliers
     * @param core
     */
    explicit SubGradient(SubInstance& subinst_) : subinst(subinst_) { }

    static LocalMultipliers u_greedy_init(const SubInstance& subinst) {

        LocalMultipliers u_0(subinst.get_nrows(), REAL_MAX);

        // init multipliers
        for (auto i = 0UL; i < subinst.get_nrows(); ++i) {
            for (const auto j : subinst.get_row(i)) {
                assert(!subinst.get_rows().empty());
                const auto candidate = subinst.get_col(j).get_cost() / static_cast<real_t>(subinst.get_col(j).size());
                if (candidate < u_0[i]) { u_0[i] = candidate; }
            }
        }

        return u_0;
    }

    static LocalMultipliers u_perturbed_init(const LocalMultipliers& u_star, std::mt19937& rnd) {

        auto u_0 = LocalMultipliers(u_star.size());
        auto dist = std::uniform_real_distribution<real_t>(0.9, 1.1);

        idx_t u0size = u_0.size();
        for (idx_t i = 0; i < u0size; i++) { u_0[i] = dist(rnd) * u_star[i]; }

        return u_0;
    }

    LocalMultipliers solve(real_t UB, const LocalMultipliers& u_0, Timer& time_limit) {
        size_t max_iter = 10 * subinst.get_rows().size();
        idx_t nrows = subinst.get_rows().size();

        // PARAMETERS
        lambda = StepSizeFactor(0.1, 20);                         // step size
        PricingPeriod T(10, std::min<idx_t>(1000UL, nrows / 3));  // pricing frequency
        ExitCondition exit_now(300U);

        LB_star = REAL_LOWEST;
        u = u_0;
        u_star = u_0;

        MStar covered_rows;
        LocalSolution S;
        real_t real_LB = lb_maintainer.compute(subinst, u);

        std::vector<std::pair<idx_t, real_t>> delta_u;

        for (idx_t iter = 0; iter < max_iter; ++iter) {
            if (time_limit.exceeded_tlim()) { break; }

            lambda.update(real_LB);
            norm_reducer.compute_reduced_sol(subinst, S, covered_rows);

            // Multipliers update:
            delta_u.clear();
            idx_t s2sum = 0;
            for (idx_t cov : covered_rows) { s2sum += (1 - static_cast<int>(cov)) * (1 - static_cast<int>(cov)); }

            if (s2sum > 0) {
                for (idx_t i = 0; i < nrows; ++i) {
                    real_t new_u = std::max<real_t>(0.0, u[i] + lambda.get() * ((UB - real_LB) / s2sum) * (1 - static_cast<int>(covered_rows[i])));
                    if (std::abs(new_u - u[i]) > REAL_TOLERANCE) {
                        delta_u.emplace_back(i, new_u - u[i]);
                        u[i] = new_u;
                    }
                }
            } else {
                IF_VERBOSE { fmt::print(" WARNING: s2sum == 0\n"); }
                return u_star;
            }

            real_LB = lb_maintainer.update(subinst, delta_u);

            if (real_LB > LB_star) {
                LB_star = real_LB;
                u_star = u;
            }

            // fmt::print("[{:^4}] Lower Bound: {:.4} (best {:.4}), lambda {:.4}, S_cost {:.4}\n", iter, real_LB, LB_star, lambda.get(),
            // S.compute_cost(subinst));

            if (real_LB > UB - HAS_INTEGRAL_COSTS) {
                IF_VERBOSE { fmt::print(" WARNING: real_LB({}) > UB({}) - {}\n", real_LB, UB, HAS_INTEGRAL_COSTS); }
                u_star = u;
                return u_star;
            }

            if (covered_rows.get_uncovered() == 0) {
                real_t S_cost = S.compute_cost(subinst);
                if (S_cost < UB) { UB = S_cost; }
            }

            if (exit_now(LB_star)) { return u_star; }

            T.inc();
            if (T.reached()) {
                const auto global_LB = subinst.price(u);
                T.reset(global_LB, real_LB, UB);
                // fmt::print("Pricing: global {}, local {}\n", global_LB, real_LB);

                LB_star = real_LB = lb_maintainer.compute(subinst, u);
            }
        }
        return u_star;
    }


    real_t get_best_LB() { return LB_star; }

private:
    SubInstance& subinst;

    real_t LB_star;
    LocalMultipliers u_star;
    LocalMultipliers u;

    ReducedLagrMultLB<1, 1000> norm_reducer;
    DeltaLowerBound<0, 1> lb_maintainer;

    StepSizeFactor lambda;
};

#endif


/* ################################################################# */
/* #### Original Header: include/TwoPhase.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_TWOPHASE_HPP_
#define SPH_INCLUDE_TWOPHASE_HPP_

#include <fmt/core.h>

#include <algorithm>
#include <cassert>
#include <vector>

/* #include "ExactSolver.hpp" */
/* #include "Multipliers.hpp" */
/* #include "Solution.hpp" */
/* #include "SubGradient.hpp" */
/* #include "SubInstance.hpp" */
/* #include "cft.hpp" */

class TwoPhase {
public:
    TwoPhase(SubInstance& subinst_) : subinst(subinst_), subgradient(subinst_) { }

    inline GlobalSolution solve(const real_t global_UB, GlobalSolution& S_star, Timer& exact_time_limit) {
        return operator()(global_UB, S_star, exact_time_limit);
    }

    GlobalSolution operator()(const real_t global_UB, GlobalSolution& S_star, Timer& exact_time_limit) {

        LocalMultipliers u_star(subinst.get_nrows(), 0.0);

        real_t glb_UB_star = std::min<real_t>(global_UB, S_star.get_cost());
        real_t fixed_cost = subinst.get_fixed_cost();
        real_t subgrad_UB = glb_UB_star - fixed_cost;

        // 1. SUBGRADIENT PHASE
        LocalMultipliers u_k = subgradient.solve(subgrad_UB, SubGradient::u_greedy_init(subinst), subinst.get_timelimit());

        real_t lcl_LB = subgradient.get_best_LB();
        real_t glb_LB = fixed_cost + lcl_LB;

        if (fixed_cost == 0.0) { glo_u = GlobalMultipliers(subinst, u_k); }
        if (glb_LB >= glb_UB_star - HAS_INTEGRAL_COSTS) { return S_star; }

        // 2. HEURISTIC PHASE
        LocalSolution S_curr(subinst.get_localized_solution(S_star));
        cplex_heur(S_curr, S_star, glb_UB_star, exact_time_limit);

        return S_star;
    }

    inline GlobalMultipliers& get_global_u() { return glo_u; }

private:
    void cplex_heur(LocalSolution& S_init, GlobalSolution& S_star, real_t& glb_UB_star, Timer& exact_time_limit) {

        real_t fixed_cost = subinst.get_fixed_cost();
        real_t S_init_cost = S_init.compute_cost(subinst);
        IF_VERBOSE { fmt::print("│ Initial solution value {} (global: {})\n", S_init_cost, S_init_cost + fixed_cost); }

        LocalSolution S = exact.build_and_opt(subinst, S_init, exact_time_limit);

        if (S.size() > 0) {

            real_t S_cost = S.compute_cost(subinst);

            real_t gS_cost = fixed_cost + S_cost;
            subinst.update_sol_costs(S, gS_cost);

            if (gS_cost < glb_UB_star) {

                glb_UB_star = gS_cost;
                S_star = GlobalSolution(subinst, S);
                IF_VERBOSE { fmt::print("│ ══> CPLEX improved global UB: {} (fixed {} + local-cost {})\n", S_star.get_cost(), fixed_cost, S_cost); }

            } else {
                IF_VERBOSE { fmt::print("│ ──> CPLEX Improved local UB: {} (global value {}, best is {})\n", S_cost, S_cost + fixed_cost, glb_UB_star); }
            }
        }

        IF_VERBOSE { fmt::print("└───────────────────────────────────────────────────────────────────────────────────\n\n"); }
    }


private:
    SubInstance& subinst;

    SubGradient subgradient;
    ExactSolver exact;

    GlobalMultipliers glo_u;
};


#endif

/* ################################################################# */
/* #### Original Header: include/Refinement.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_REFINEMENT_HPP_
#define SPH_INCLUDE_REFINEMENT_HPP_

/* #include "Instance.hpp" */
/* #include "SubInstance.hpp" */
/* #include "TwoPhase.hpp" */
/* #include "cft.hpp" */

#define ALPHA 1.1
#define BETA 0.8
#define PI_MIN 0.3

#define PI_MAX 0.9
#define POST_OPT_TRIALS 100

#define SHORT_T_LIM 10.0
#define LONG_T_LIM(TOTAL_TIME) (std::min(TOTAL_TIME / 2.0, 100.0))

class Refinement {
public:
    Refinement(Instance& inst_, std::mt19937& rnd_) : inst(inst_), subinst(inst_), two_phase(subinst), rnd(rnd_) { }

    // S_init must be a global complete solution
    template <unsigned long ROUTES_HARD_CAP>
    std::vector<idx_t> solve([[maybe_unused]] const std::vector<idx_t>& S_init) {

        // 1.
        GlobalSolution S_star;

        if (!S_init.empty()) {
            real_t cost = 0.0;
            for (idx_t j : S_init) { cost += inst.get_col(j).get_cost(); }
            for (auto j : S_init) { S_star.push_back(j); }
            S_star.set_cost(cost);
            IF_VERBOSE { fmt::print("Found warm start with cost {}.\n", cost); }
        }

        GlobalMultipliers u_star;

        real_t pi = PI_MIN;
        real_t last_improving_pi = pi;

        int post_optimization_trials = POST_OPT_TRIALS;
        Timer& global_time_limit = inst.get_timelimit();

        idx_t iter = 1;
        do {
            subinst.reset();  // 2.

            {
                Timer iteration_timer = Timer(iter == 1 ? LONG_T_LIM(global_time_limit.seconds_until_end()) : SHORT_T_LIM);
                GlobalSolution S = two_phase(S_star.get_cost(), S_star, iteration_timer);  // 3. & 4.

                assert(!(std::fabs(pi - PI_MIN) > 0.001 && inst.get_fixed_cols().empty()));
                pi *= ALPHA;  // 6.

                if (S.get_cost() < S_star.get_cost()) {  // update best solution
                    S_star = std::move(S);               // 5.

                    last_improving_pi = pi;

                    pi = std::max(pi / (ALPHA * ALPHA), PI_MIN);  // 6.
                }

                if (S_star.get_cost() - 1.0 <= BETA * u_star.get_lb() || pi > PI_MAX || inst.get_active_rows_size() <= 0) {

                    if (post_optimization_trials <= 0) {
                        IF_VERBOSE {
                            fmt::print("╔═ REFINEMENT: iter {:2} ═════════════════════════════════════════════════════════════\n", iter);
                            fmt::print("║ Early Exit: β(={:.1f}) * LB(={:.1f}) > UB(={:.1f}) - 1\n", BETA, u_star.get_lb(), S_star.get_cost());
                            fmt::print("║ Active rows {}, active cols {}, pi {:.3f}\n", inst.get_active_rows_size(), inst.get_active_cols().size(), pi);
                            fmt::print("║ LB {:.1f}, UB {:.1f}, UB size {}\n", u_star.get_lb(), S_star.get_cost(), S_star.size());
                            fmt::print("╚═══════════════════════════════════════════════════════════════════════════════════\n\n");
                        }
                        break;
                    }

                    --post_optimization_trials;
                    fmt::print("   POST-OPTIMIZATION REFINEMENT: iter {:2}\n", POST_OPT_TRIALS - post_optimization_trials);

                    pi = last_improving_pi;
                    last_improving_pi = std::max(PI_MIN, last_improving_pi / ALPHA);
                }
            }  //(S and u are potentially been moved, better to encapsulate them into a block)

            if (iter == 1) {
                u_star = two_phase.get_global_u();
                std::vector<idx_t> old_to_new_idx_map = inst.prune_instance<ROUTES_HARD_CAP>(u_star);

                IF_VERBOSE { fmt::print("Instance size: {}x{}\n", inst.get_nrows(), inst.get_ncols()); }

                if (!old_to_new_idx_map.empty()) {
                    for (idx_t& gj : S_star) {
                        gj = old_to_new_idx_map[gj];
                        assert(gj != REMOVED_INDEX);
                    }
                }
            }

            inst.reset_fixing();

            if (global_time_limit.exceeded_tlim()) {
                IF_VERBOSE {
                    fmt::print("╔═ REFINEMENT: iter {:2} ═════════════════════════════════════════════════════════════\n", iter);
                    fmt::print("║ Timelimit exceeded\n");
                    fmt::print("╚═══════════════════════════════════════════════════════════════════════════════════\n\n");
                }
                break;
            }

            pi = std::min<real_t>(PI_MAX, pi);

            // 7. Refinement Fix
            cols_to_fix = random_fix(S_star, pi);

            assert(cols_to_fix.size() <= S_star.size());

            std::sort(cols_to_fix.begin(), cols_to_fix.end());
            inst.fix_columns(cols_to_fix, covered_rows);

            IF_VERBOSE {
                fmt::print("╔═ REFINEMENT: iter {:2} ═════════════════════════════════════════════════════════════\n", iter);
                fmt::print("║ Active rows {}, active cols {}, pi {:.3f}\n", inst.get_active_rows_size(), inst.get_active_cols().size(), pi);
                fmt::print("║ LB {:.1f}, UB {:.1f}, UB size {}\n", u_star.get_lb(), S_star.get_cost(), S_star.size());
                fmt::print("╚═══════════════════════════════════════════════════════════════════════════════════\n\n");
            }

            assert(inst.get_fixed_cost() <= S_star.get_cost());
            ++iter;
        } while (true);

        return S_star;
    }

private:
    [[nodiscard]] std::vector<idx_t> refinement_fix(GlobalSolution S_star, GlobalMultipliers u_star, real_t pi) {
        UniqueColSet& cols = inst.get_cols();

        covered_rows.reset_covered(cols, S_star, inst.get_nrows());

        deltas.resize(S_star.size());
        for (idx_t j = 0; j < S_star.size(); ++j) {
            auto& col = cols[S_star[j]];
            deltas[j].first = S_star[j];
            deltas[j].second = std::max<real_t>(col.compute_lagr_cost(u_star), 0.0);
            for (auto i : col) { deltas[j].second += u_star[i] * (covered_rows[i] - 1.0) / covered_rows[i]; }
        }
        std::sort(deltas.begin(), deltas.end(), [](auto& a, auto& b) { return a.second < b.second; });

        covered_rows.reset_uncovered(inst.get_nrows());
        real_t covered_fraction = 0.0;
        idx_t n = 0;
        for (; n < deltas.size() && covered_fraction < pi; ++n) {
            covered_rows.cover_rows(cols[deltas[n].first]);
            covered_fraction = static_cast<real_t>(covered_rows.get_covered()) / static_cast<real_t>(inst.get_nrows());
        }

        cols_to_fix.resize(n);
        for (idx_t j2 = 0; j2 < n; ++j2) { cols_to_fix[j2] = deltas[j2].first; }

        return cols_to_fix;
    }

    [[nodiscard]] std::vector<idx_t> random_fix(GlobalSolution S_star, real_t pi) {
        UniqueColSet& cols = inst.get_cols();

        std::shuffle(S_star.begin(), S_star.end(), rnd);

        cols_to_fix.clear();
        covered_rows.reset_uncovered(inst.get_nrows());
        real_t covered_fraction = 0.0;

        for (idx_t n = 0; n < S_star.size() && covered_fraction < pi; ++n) {
            cols_to_fix.emplace_back(S_star[n]);
            covered_rows.cover_rows(cols[S_star[n]]);
            covered_fraction = static_cast<real_t>(covered_rows.get_covered()) / static_cast<real_t>(inst.get_nrows());
        }

        return cols_to_fix;
    }

    [[nodiscard]] std::vector<idx_t> binary_tournament_fix(GlobalSolution S_star, GlobalMultipliers u_star, real_t pi) {
        UniqueColSet& cols = inst.get_cols();

        covered_rows.reset_covered(cols, S_star, inst.get_nrows());

        deltas.resize(S_star.size());
        for (idx_t j = 0; j < S_star.size(); ++j) {
            auto& col = cols[S_star[j]];
            deltas[j].first = S_star[j];
            deltas[j].second = std::max<real_t>(col.compute_lagr_cost(u_star), 0.0);
            for (auto i : col) { deltas[j].second += u_star[i] * (covered_rows[i] - 1.0) / covered_rows[i]; }
        }

        idx_t dsize = deltas.size();

        cols_to_fix.resize(dsize);
        covered_rows.reset_uncovered(inst.get_nrows());
        real_t covered_fraction = 0.0;
        idx_t n = 0;
        for (; n < deltas.size() && dsize > 1 && covered_fraction < pi; ++n, --dsize) {
            auto& cand1 = deltas[rnd() % dsize];
            idx_t c2;
            do { c2 = rnd() % dsize; } while (cand1.first == deltas[c2].first);
            auto& winner = cand1.second < deltas[c2].second ? cand1 : deltas[c2];

            cols_to_fix[n] = winner.first;
            covered_rows.cover_rows(cols[winner.first]);
            covered_fraction = static_cast<real_t>(covered_rows.get_covered()) / static_cast<real_t>(inst.get_nrows());

            std::swap(deltas.back(), winner);
        }
        cols_to_fix.resize(n);

        return cols_to_fix;
    }

    [[nodiscard]] std::vector<idx_t> random_fix2(GlobalSolution S_star, GlobalMultipliers u_star, real_t pi) {
        UniqueColSet& cols = inst.get_cols();
        std::uniform_real_distribution<real_t> dist(0.5, 1.5);

        covered_rows.reset_covered(cols, S_star, inst.get_nrows());

        deltas.resize(S_star.size());
        for (idx_t j = 0; j < S_star.size(); ++j) {
            auto& col = cols[S_star[j]];
            deltas[j].first = S_star[j];
            deltas[j].second = std::max<real_t>(col.compute_lagr_cost(u_star), 0.0);
            for (auto i : col) { deltas[j].second += u_star[i] * (covered_rows[i] - 1.0) / covered_rows[i]; }
            deltas[j].second *= dist(rnd);
        }
        std::sort(deltas.begin(), deltas.end(), [](auto& a, auto& b) { return a.second < b.second; });

        covered_rows.reset_uncovered(inst.get_nrows());
        real_t covered_fraction = 0.0;
        idx_t n = 0;
        for (; n < deltas.size() && covered_fraction < pi; ++n) {
            covered_rows.cover_rows(cols[deltas[n].first]);
            covered_fraction = static_cast<real_t>(covered_rows.get_covered()) / static_cast<real_t>(inst.get_nrows());
        }

        cols_to_fix.resize(n);
        for (idx_t j2 = 0; j2 < n; ++j2) { cols_to_fix[j2] = deltas[j2].first; }

        return cols_to_fix;
    }

private:
    Instance& inst;

    SubInstance subinst;
    MStar covered_rows;
    TwoPhase two_phase;

    std::mt19937& rnd;

    // retain allocated memory (anti-RAII, should be local)
    std::vector<std::pair<idx_t, real_t>> deltas;
    std::vector<idx_t> cols_to_fix;
};

#endif

/* ################################################################# */
/* #### Original Header: include/SPHeuristic.hpp */
/* ################################################################# */

#ifndef SPH_INCLUDE_SPHEURISTIC_HPP_
#define SPH_INCLUDE_SPHEURISTIC_HPP_

#include <cassert>

/* #include "Instance.hpp" */
/* #include "Refinement.hpp" */

class SPHeuristic {
public:
    explicit SPHeuristic(const idx_t nrows_) : inst(nrows_), rnd(std::mt19937()), refinement(inst, rnd) { }
    SPHeuristic(const idx_t nrows_, int seed) : inst(nrows_), rnd(std::mt19937(seed)), refinement(inst, rnd) { }

    [[nodiscard]] inline idx_t get_ncols() const { return inst.get_ncols(); }
    [[nodiscard]] inline idx_t get_nrows() const { return inst.get_nrows(); }
    [[nodiscard]] inline UniqueColSet &get_cols() { return inst.get_cols(); }
    [[nodiscard]] inline Column &get_col(idx_t idx) { return inst.get_col(idx); }
    [[nodiscard]] inline const Column &get_col(idx_t idx) const { return inst.get_col(idx); }

    void inline set_timelimit(double seconds) { inst.set_timelimit(seconds); }

    [[nodiscard]] inline Timer &get_timelimit() { return inst.get_timelimit(); }

    template <typename UniqueColContainer>
    std::vector<idx_t> inline add_columns(const UniqueColContainer &new_cols) {
        return inst.add_columns(new_cols);
    }

    std::vector<idx_t> inline add_columns(const std::vector<real_t> &costs, const std::vector<real_t> &sol_costs, const std::vector<idx_t> &matbeg,
                                          const std::vector<idx_t> &matval) {
        return inst.add_columns(costs, sol_costs, matbeg, matval);
    }

    template <unsigned long ROUTES_HARD_CAP = INST_HARD_CAP>
    std::vector<idx_t> inline solve([[maybe_unused]] const std::vector<idx_t> &S_init) {
        return refinement.solve<ROUTES_HARD_CAP>(S_init);
    }

private:
    Instance inst;
    std::mt19937 rnd;
    Refinement refinement;
};

#endif
