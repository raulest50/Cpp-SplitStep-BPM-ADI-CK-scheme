
#include "LaserBeam.h"


std::vector< std::complex<float> >
LaserBeam::GaussianBeamSaleh2D(
        std::size_t nx, std::size_t ny,
        float dx, float dy,
        float waist, float wavelength, float z, float E0
        ) {

        const float pi = 3.14159265358979323846f;
        float k  = 2.0f * pi / wavelength;                // propagation constant
        float z0 = pi * waist * waist / wavelength;      // Rayleigh range

        // beam radius at z
        float wz = waist * std::sqrt(1.0f + (z*z)/(z0*z0));
        // radius of curvature at z, avoid division by zero
        float Rz = (z == 0.0f)
                   ? std::numeric_limits<float>::infinity()
                   : z * (1.0f + (z0*z0)/(z*z));
        float psi = std::atan(z / z0);                   // Gouy phase

        std::vector<std::complex<float>> beam;
        beam.reserve(nx * ny);

        // center the grid around zero
        float x0 = (static_cast<float>(nx) - 1.0f)/2.0f;
        float y0 = (static_cast<float>(ny) - 1.0f)/2.0f;

        for (std::size_t i = 0; i < nx; ++i) {
                float x = (static_cast<float>(i) - x0) * dx;
                for (std::size_t j = 0; j < ny; ++j) {
                        float y = (static_cast<float>(j) - y0) * dy;
                        float r2 = x*x + y*y;

                        // amplitude term
                        float amp = E0 * (waist / wz)
                                    * std::exp(-r2/(wz*wz));

                        // total phase
                        float phase = -k*z
                                      - (k * r2)/(2.0f * Rz)
                                      + psi;

                        // polar form gives amp·e^{i·phase}
                        beam.emplace_back(std::polar(amp, phase));
                }
        }
        return beam;
}
