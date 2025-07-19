
#include "Chemkin/reactionConstants.hh"

#include <cmath>

#include "Chemkin/constants.hh"

auto lab109::chemkin::reaction_constants::areniusConstant(double A,
                                                          double N,
                                                          double E,
                                                          double T) -> double {
    return A * std::pow(T, N) * exp(-E / T);
}

auto lab109::chemkin::reaction_constants::equilibriumConstant(std::span<const int>    a,
                                                              std::span<const int>    b,
                                                              std::span<const double> H,
                                                              std::span<const double> S,
                                                              double                  T) -> double {
    int    sum_b_m_a = {};
    double dH        = {};
    double dS        = {};

    for (std::size_t i = {}; i < a.size(); ++i) {
        const int b_m_a = b[i] - a[i];
        dH += b_m_a * H[i];
        dS += b_m_a * S[i];
        sum_b_m_a += b_m_a;
    }

    const auto k_p = std::exp((dS - dH / T) / constants::K_R_GAS);
    return k_p * std::pow(constants::K_PA_TO_ATM / constants::K_R_GAS / T, sum_b_m_a);
}