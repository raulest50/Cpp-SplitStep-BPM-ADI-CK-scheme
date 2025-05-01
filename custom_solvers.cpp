//
// Created by Usuario on 01/05/2025.
//

#include "custom_solvers.h"

void custom_solvers::thomas1(
        const std::vector<std::complex<float>>& a,  // Sub-diagonal (a[0] unused)
        const std::vector<std::complex<float>>& b,  // Main diagonal
        const std::vector<std::complex<float>>& c,  // Super-diagonal (c[N-1] unused)
        const std::vector<std::complex<float>>& d,  // Right-hand side
        std::vector<std::complex<float>>& x         // Solution
) {
    int N = b.size();
    std::vector<std::complex<float>> c_prime(N);
    std::vector<std::complex<float>> d_prime(N);

    // First modified coefficient
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    // Forward sweep
    for (int i = 1; i < N; ++i) {
        std::complex<float> denom = b[i] - a[i] * c_prime[i - 1];
        c_prime[i] = (i < N - 1) ? c[i] / denom : 0.0f;
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom;
    }

    // Back substitution
    x[N - 1] = d_prime[N - 1];
    for (int i = N - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }
}