//
// Created by Usuario on 01/05/2025.
//

#ifndef CPPSPLITSTEPBPMADICKSCHEME_CUSTOM_SOLVERS_H
#define CPPSPLITSTEPBPMADICKSCHEME_CUSTOM_SOLVERS_H

#include <vector>
#include <complex>

class custom_solvers {
public:
    // Nueva función personalizada con estructura especial
    static void custom_thomas_solver(
            const std::complex<float>& dp,   // Valor para todos los elementos en la diagonal principal excepto primero y último
            const std::complex<float>& dp1,  // Valor para el primer elemento en la diagonal principal [0,0]
            const std::complex<float>& dp2,  // Valor para el último elemento en la diagonal principal [-1,-1]
            const std::complex<float>& do_,  // Valor para todos los elementos en las diagonales secundarias
            const std::vector<std::complex<float>>& b,  // Vector del lado derecho
            std::vector<std::complex<float>>& x         // Vector solución
    );
    
    // Multiplica una matriz tridiagonal con estructura especial por un vector
    static std::vector<std::complex<float>> compute_b_vector(
            const std::complex<float>& dp,   // Valor para todos los elementos en la diagonal principal excepto primero y último
            const std::complex<float>& dp1,  // Valor para el primer elemento en la diagonal principal [0,0]
            const std::complex<float>& dp2,  // Valor para el último elemento en la diagonal principal [-1,-1]
            const std::complex<float>& do_,  // Valor para todos los elementos en las diagonales secundarias
            const std::vector<std::complex<float>>& x0   // Vector de entrada a multiplicar
    );
    
    // Realiza un paso ADI en la dirección x
    static std::vector<std::complex<float>> adi_x(
            const std::vector<std::complex<float>>& phi,  // Campo complejo 2D aplanado (row-major)
            int Nx,                                       // Número de puntos en dirección x
            int Ny,                                       // Número de puntos en dirección y
            float eps,                                    // Valor pequeño para evitar división por cero
            float k,                                      // Número de onda
            float dz,                                     // Paso en dirección z
            float dx                                      // Paso en dirección x
    );
    
    // Realiza un paso ADI en la dirección y
    static std::vector<std::complex<float>> adi_y(
            const std::vector<std::complex<float>>& phi,  // Campo complejo 2D aplanado (row-major)
            int Nx,                                       // Número de puntos en dirección x
            int Ny,                                       // Número de puntos en dirección y
            float eps,                                    // Valor pequeño para evitar división por cero
            float k,                                      // Número de onda
            float dz,                                     // Paso en dirección z
            float dy                                      // Paso en dirección y
    );
    
    // Aplica el operador no lineal (medio paso)
    static std::vector<std::complex<float>> half_nonlinear(
            const std::vector<std::complex<float>>& phi,  // Campo complejo
            float k_sample,                               // Número de onda en el medio
            float n2_sample,                              // Índice de refracción no lineal
            float dz                                      // Paso en dirección z
    );
    
    // Aplica la absorción lineal (medio paso)
    static std::vector<std::complex<float>> half_linear_absorption(
            const std::vector<std::complex<float>>& phi,  // Campo complejo
            float alpha,                                  // Coeficiente de absorción lineal
            float dz                                      // Paso en dirección z
    );
    
    // Aplica la absorción de dos fotones (medio paso)
    static std::vector<std::complex<float>> half_2photon_absorption(
            const std::vector<std::complex<float>>& phi,  // Campo complejo
            float beta,                                   // Coeficiente de absorción de dos fotones
            float dz                                      // Paso en dirección z
    );
};

#endif //CPPSPLITSTEPBPMADICKSCHEME_CUSTOM_SOLVERS_H