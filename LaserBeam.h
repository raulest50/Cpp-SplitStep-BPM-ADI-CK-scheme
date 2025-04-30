
#ifndef LASERBEAM_H  // include guard start
#define LASERBEAM_H
#include <complex>
#include <vector>


class LaserBeam {
public:

    static std::vector< std::complex<float> >
    GaussianBeamSaleh2D(
        std::size_t nx,
        std::size_t ny,
        float dx, float dy,
        float waist, float wavelength,float z, float E0 = 1.0f
        );

};



#endif //LASERBEAM_H
