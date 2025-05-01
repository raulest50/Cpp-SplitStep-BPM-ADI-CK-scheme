//
// Created by Usuario on 01/05/2025.
//

#ifndef CPPSPLITSTEPBPMADICKSCHEME_CUSTOM_SOLVERS_H
#define CPPSPLITSTEPBPMADICKSCHEME_CUSTOM_SOLVERS_H

#include <vector>
#include <complex>

class custom_solvers {
public:
    static void thomas1(
            const std::vector<std::complex<float>>& a,  // Sub-diagonal (a[0] unused)
            const std::vector<std::complex<float>>& b,  // Main diagonal
            const std::vector<std::complex<float>>& c,  // Super-diagonal (c[N-1] unused)
            const std::vector<std::complex<float>>& d,  // Right-hand side
            std::vector<std::complex<float>>& x         // Solution
    );



};


#endif //CPPSPLITSTEPBPMADICKSCHEME_CUSTOM_SOLVERS_H
