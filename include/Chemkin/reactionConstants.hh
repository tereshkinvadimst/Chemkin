#ifndef LAB109_CHEMKIN_REACTION_CONSTANTS_H
#define LAB109_CHEMKIN_REACTION_CONSTANTS_H
#include <span>


namespace lab109::chemkin::reaction_constants {

auto areniusConstant(double A, double N, double E, double T) -> double;

auto equilibriumConstant(std::span<const int>    a,
                         std::span<const int>    b,
                         std::span<const double> H,
                         std::span<const double> S,
                         double                  T) -> double;

}  // namespace lab109::chemkin::reaction_constants

#endif  // LAB109_CHEMKIN_REACTION_CONSTANTS_H