#include <cstdio>
#include <complex>
#include <vector>
#include "LaserBeam.h"
#include "PrintHelpers.h"


int main() {

    size_t nx =47, ny = 47;
    float dx = 15e-6, dy = 15e-6;


    float n2 = 2.5e-20;     // second order sample refraction index
    float k;                // vacuum wave number

    float waist = 50e-5f;   // w₀ = 50 μm
    float wavelength = 800e-9;     // λ = 800 nm
    float z = 0;            // at the waist plane
    float E0 = 1.0f;        // normalized peak amplitude

    std::vector< std::complex<float> > Phi_k0 =
        LaserBeam::GaussianBeamSaleh2D(nx, ny, dx, dy, waist,wavelength, z, E0);

    PrintHelpers::printBeam2D(Phi_k0, nx, ny);

    return 0;

}






