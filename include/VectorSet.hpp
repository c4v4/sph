// Copyright (c) 2022 Francesco Cavaliere
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef CAV_VECTORSET_HPP
#define CAV_VECTORSET_HPP

#include <cstddef>
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
         * @brief Since the set does not contain T, but T*,
         * this hacky struct is used to make the T* behave like
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

            T*& base_addr;
            size_t idx;
        };

        using Set = std::unordered_set<TPtrWrap, Hash>;
        using siter = typename Set::iterator;
        using const_siter = typename Set::const_iterator;

    public:
        VectorSet() : vec_data(new T*) { }

        template <typename Iter>
        VectorSet(Iter beg, Iter end) {
            reserve(end - beg);
            for (auto it = beg; it != end; ++it) { emplace_back(*it); }
        }

        VectorSet(const VectorSet& other) : vec(other.vec), vec_data(new T*), set(other.set) { *vec_data = vec.data(); }
        VectorSet(VectorSet&& other) noexcept : vec(std::move(other.vec)), vec_data(std::exchange(other.vec_data, nullptr)), set(std::move(other.set)) { }
        ~VectorSet() { delete vec_data; }

        VectorSet& operator=(const VectorSet& other) {
            vec = other.vec;
            set = other.set;
            *vec_data = vec.data();
            return *this;
        }

        VectorSet& operator=(VectorSet&& other) noexcept {
            vec = std::move(other.vec);
            set = std::move(other.set);
            vec_data = std::exchange(other.vec_data, nullptr);
            return *this;
        }

        template <typename... _Args>
        std::pair<viter, bool> emplace_back(_Args&&... args) {
            T elem(std::forward<_Args>(args)...);
            T* elem_ptr = std::addressof(elem);
            auto wrapped_elem = TPtrWrap(elem_ptr, 0);

            auto set_elem_iter = set.find(wrapped_elem);
            if (set_elem_iter == set.end()) {
                vec.emplace_back(std::move(elem));
                *vec_data = vec.data();
                size_t elem_position = vec.size() - 1;
                set.emplace(*vec_data, elem_position);
                return std::make_pair<viter, bool>(vec.begin() + elem_position, true);
            }

            ptrdiff_t elem_position = set_elem_iter->data() - vec.data();
            return std::make_pair<viter, bool>(vec.begin() + elem_position, false);
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
            TPtrWrap wrapped_elem(elem_ptr, 0);

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

        bool is_corrupted() {
            for (auto elem : set) {
                if (elem.base_addr != *vec_data) {
                    fmt::print(stderr, "VectorSet corruped!");
                    return true;
                }
                if (elem.data() < vec.data() || elem.data() >= vec.data() + vec.size()) {
                    fmt::print(stderr, "VectorSet corruped!");
                    return true;
                }
            }
            return false;
        }

    private:
        std::vector<T> vec;
        T** vec_data = nullptr;  // we can't get vec.data() ptr ref, so this is a reproduction that needs to be kept updated

        Set set;
    };
}  // namespace cav

#endif