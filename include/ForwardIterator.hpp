#ifndef CAV_FORWARDITERATOR_HPP
#define CAV_FORWARDITERATOR_HPP

#include <iterator>
#include <numeric>

namespace cav {
    template <typename Base, class OperInc, class OperDefer>
    struct ForwardIterator {
        typedef typename std::forward_iterator_tag iterator_category;
        typedef typename std::ptrdiff_t difference_type;
        typedef typename std::remove_reference_t<typename std::result_of_t<OperDefer(Base&)>> value_type;
        typedef value_type* pointer;
        typedef value_type& reference;

    public:
        ForwardIterator(Base base_) : base(base_){};

        inline reference operator*() { return OperDefer()(base); }

        inline ForwardIterator operator+(difference_type i) noexcept {
            Base _base = base;
            for (; i > 0; --i)
                _base = OperInc()(_base);
            return ForwardIterator(_base);
        }

        inline ForwardIterator operator++() noexcept {
            base = OperInc()(base);
            return *this;
        }

        inline ForwardIterator operator+=(difference_type i) noexcept {
            for (; i > 0; --i)
                base = OperInc()(base);
            return *this;
        }

        inline ForwardIterator operator++(int) noexcept {
            ForwardIterator curr = base;
            base = OperInc()(base);
            return curr;
        }

        inline Base& data() { return base; }

        inline bool operator==(const ForwardIterator& other) const noexcept { return base == other.base; }
        inline bool operator!=(const ForwardIterator& other) const noexcept { return base != other.base; }

    private:
        Base base;
    };
}  // namespace cav

#endif