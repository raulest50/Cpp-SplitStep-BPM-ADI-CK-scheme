#ifndef LASERBEAM_H  // include guard start
#define LASERBEAM_H
#include <complex>
#include <vector>


class LaserBeam {
public:
    // Constantes del láser
    static constexpr float WAVELENGTH = 800e-9f;  // m - longitud de onda (800 nm)
    static constexpr float WAIST = 3e-6f;         // m - beam waist (3 μm)
    static constexpr float IPEAK = 1e10f;        // W/m² - intensidad pico
    
    /**
     * Genera un campo TEM00 gaussiano complejo E(x,y) listo para usar en BPM.
     * 
     * @param nx Número de puntos en dirección x
     * @param ny Número de puntos en dirección y
     * @param dx Paso espacial en dirección x (metros)
     * @param dy Paso espacial en dirección y (metros)
     * @param w0 Radio del waist (haz) en metros
     * @param I0 Intensidad pico en W/m²
     * @param fase_inicial Fase global opcional (en radianes)
     * @return Vector con campo eléctrico complejo E(x,y) aplanado (row-major)
     */
    static std::vector<std::complex<float>> campo_tem00(
            size_t nx, size_t ny,
            float dx, float dy,
            float w0, 
            float I0,
    );
};

#endif //LASERBEAM_H