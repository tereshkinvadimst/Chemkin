#ifndef LAB109_KINETICS_HH
#define LAB109_KINETICS_HH
#include <cstdint>
#include <mdspan>
#include <span>

namespace lab109::chemkin {

void computeKineticsParameters(std::mdspan<const std::int8_t, std::dextents<std::size_t, 2>> a,
                               std::mdspan<const std::int8_t, std::dextents<std::size_t, 2>> b,
                               std::span<const double>                                       k_f,
                               std::span<const double>                                       k_r,
                               std::span<const double>                                       X,
                               std::span<double>                                             dW_f,
                               std::span<double>                                             dW_r,
                               std::span<double>                                             dX_f,
                               std::span<double>                                             dX_r);

}  // namespace lab109::chemkin

#endif  // LAB109_KINETICS_HH