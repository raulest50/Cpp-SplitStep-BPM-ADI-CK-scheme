#include "LaserBeam.h"
#include <cmath>

std::vector<std::complex<float>> LaserBeam::campo_tem00(
        size_t nx, size_t ny,
        float dx, float dy,
        float w0, 
        float I0,
) {
    // Crear el vector para el campo resultante
    std::vector<std::complex<float>> campo(nx * ny);
    
    // Calcular el campo en cada punto de la malla
    for (size_t i = 0; i < nx; ++i) {
        // Coordenada x (centrada)
        float x = (static_cast<float>(i) - static_cast<float>(nx) / 2.0f) * dx;
        
        for (size_t j = 0; j < ny; ++j) {
            // Coordenada y (centrada)
            float y = (static_cast<float>(j) - static_cast<float>(ny) / 2.0f) * dy;
            
            // Calcular R² = x² + y²
            float R2 = x * x + y * y;
            
            // Amplitud gaussiana
            float Ex = std::sqrt(I0) * std::exp(-R2 / (w0 * w0));
            
            // Aplicar factor de fase
            campo[i * ny + j] = Ex;
        }
    }
    
    return campo;
}