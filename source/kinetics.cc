#include "Chemkin/kinetics.hh"

#include <cmath>

void lab109::chemkin::computeKineticsParameters(
    std::mdspan<const std::int8_t, std::dextents<std::size_t, 2>> a,
    std::mdspan<const std::int8_t, std::dextents<std::size_t, 2>> b,
    std::span<const double>                                       k_f,
    std::span<const double>                                       k_r,
    std::span<const double>                                       X,
    std::span<double>                                             dW_f,
    std::span<double>                                             dW_r,
    std::span<double>                                             dX_f,
    std::span<double>                                             dX_r) {
    const auto n_species   = X.size();
    const auto n_reactions = dW_f.size();

    for (std::size_t i = {}; i < n_reactions; ++i) {
        dW_f[i] = k_f[i];
        dW_r[i] = k_r[i];

        for (std::size_t k = {}; k < n_species; ++k) {
            dW_f[i] *= std::pow(X[k], a[i, k]);
            dW_r[i] *= std::pow(X[k], b[i, k]);
        }
    }

    for (std::size_t k = {}; k < n_species; ++k) {
        dX_f[k] = 0.;
        dX_r[k] = 0.;

        for (std::size_t i = {}; i < n_reactions; ++i) {
            const auto b_m_a = b[i, k] - a[i, k];
            dX_f[k] += b_m_a * dW_f[i];
            dX_r[k] += b_m_a * dW_r[i];
        }
    }
}