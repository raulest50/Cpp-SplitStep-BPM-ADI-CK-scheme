

#include "PrintHelpers.h"
#include <iomanip>
#include <iostream>

// Prints a row-major flat vector as an nx×ny grid
void PrintHelpers::printBeam2D(const std::vector<std::complex<float>>& beam,
                 std::size_t nx, std::size_t ny)
{
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            const auto& c = beam[i*ny + j];
            // default format for complex<T> is "(real,imag)"
            // you can also use c.real() and c.imag() if you want custom formatting:
            std::cout
                << std::fixed
                << std::setprecision(2)
                << std::setw(8) << c.real()
                << (c.imag() < 0 ? " - " : " + ")
                << std::setw(7) << std::abs(c.imag()) << "i  ";
        }
        std::cout << "\n";
    }
}