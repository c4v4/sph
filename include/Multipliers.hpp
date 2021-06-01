#ifndef AC_CFT_INCLUDE_MULTIPLIERS_HPP_
#define AC_CFT_INCLUDE_MULTIPLIERS_HPP_
#include <vector>

#include "Instance.hpp"
#include "LowerBound.hpp"
#include "cft.hpp"

class LocalMultipliers;
class GlobalMultipliers;

class BaseMultipliers : public std::vector<real_t> {
public:
    template <typename... Args>
    explicit BaseMultipliers(Args&&... _args) : std::vector<real_t>(std::forward<Args>(_args)...) { }

protected:
    template <typename Inst>
    [[nodiscard]] inline real_t compute_lb(Inst& inst) const {
        auto& cols = inst.get_cols();
        real_t LB = std::reduce(begin(), end(), static_cast<real_t>(0.0));
        for (auto& col : cols) {
            real_t c_u = col.compute_lagr_cost(*this);
            if (c_u < 0.001) { LB += c_u; }
        }

        return LB;
    }
};


class LocalMultipliers : public BaseMultipliers {
public:
    template <typename... Args>
    explicit LocalMultipliers(Args&&... _args) : BaseMultipliers(std::forward<Args>(_args)...) { }

    [[nodiscard]] inline real_t compute_lb(SubInstance& subinst) const { return BaseMultipliers::compute_lb(subinst); }
};


class GlobalMultipliers : public BaseMultipliers {
public:
    GlobalMultipliers() : lb(std::numeric_limits<real_t>::min()) { }
    explicit GlobalMultipliers(idx_t size) : BaseMultipliers(size), lb(std::numeric_limits<real_t>::min()) { }
    GlobalMultipliers(idx_t size, real_t val) : BaseMultipliers(size, val), lb(std::numeric_limits<real_t>::min()) { }
    GlobalMultipliers(const GlobalMultipliers& other) : BaseMultipliers(other), lb(other.lb) { }
    GlobalMultipliers(GlobalMultipliers&& other)  noexcept : BaseMultipliers(std::move(other)), lb(other.lb) { }

    GlobalMultipliers& operator=(const GlobalMultipliers& other) = default;
    GlobalMultipliers& operator=(GlobalMultipliers&& other) = default;

    GlobalMultipliers(SubInstance& subinst, LocalMultipliers mult) : BaseMultipliers(subinst.get_instance().get_nrows(), static_cast<real_t>(0.0)) {

        for (idx_t i = 0; i < mult.size(); ++i) {
            (*this)[subinst.get_global_row_idx(i)] = mult[i];
        }
        lb = compute_lb(subinst.get_instance());
    }

    [[nodiscard]] inline real_t get_lb() const { return lb; }

private:
    real_t lb;
};


#endif