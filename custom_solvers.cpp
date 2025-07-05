//
// Created by Usuario on 01/05/2025.
//

#include "custom_solvers.h"

// Implementación del solucionador de Thomas con estructura especial
void custom_solvers::custom_thomas_solver(
        const std::complex<float>& dp,   // Valor para todos los elementos en la diagonal principal excepto primero y último
        const std::complex<float>& dp1,  // Valor para el primer elemento en la diagonal principal [0,0]
        const std::complex<float>& dp2,  // Valor para el último elemento en la diagonal principal [-1,-1]
        const std::complex<float>& do_,  // Valor para todos los elementos en las diagonales secundarias
        const std::vector<std::complex<float>>& b,  // Vector del lado derecho
        std::vector<std::complex<float>>& x        // Vector solución
) {
    int n = b.size();

    // Crear arreglos para los coeficientes modificados
    std::vector<std::complex<float>> c_prime(n - 1);  // Diagonal superior
    std::vector<std::complex<float>> d_prime(n);      // Lado derecho modificado

    // Eliminación hacia adelante
    // Primera fila
    c_prime[0] = do_ / dp1;
    d_prime[0] = b[0] / dp1;

    // Filas intermedias
    for (int i = 1; i < n - 1; ++i) {
        std::complex<float> denominator = dp - do_ * c_prime[i - 1];
        c_prime[i] = do_ / denominator;
        d_prime[i] = (b[i] - do_ * d_prime[i - 1]) / denominator;
    }

    // Última fila
    d_prime[n - 1] = (b[n - 1] - do_ * d_prime[n - 2]) / (dp2 - do_ * c_prime[n - 2]);

    // Sustitución hacia atrás
    x[n - 1] = d_prime[n - 1];

    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }
}


// Multiplica una matriz tridiagonal con estructura especial por un vector
std::vector<std::complex<float>> custom_solvers::compute_b_vector(
        const std::complex<float>& dp,   // Valor para todos los elementos en la diagonal principal excepto primero y último
        const std::complex<float>& dp1,  // Valor para el primer elemento en la diagonal principal [0,0]
        const std::complex<float>& dp2,  // Valor para el último elemento en la diagonal principal [-1,-1]
        const std::complex<float>& do_,  // Valor para todos los elementos en las diagonales secundarias
        const std::vector<std::complex<float>>& x0   // Vector de entrada a multiplicar
) {
    int n = x0.size();
    std::vector<std::complex<float>> b(n);
    
    // Primera fila: b[0] = dp1 * x0[0] + do * x0[1]
    b[0] = dp1 * x0[0] + do_ * x0[1];
    
    // Filas intermedias: b[i] = do * x0[i-1] + dp * x0[i] + do * x0[i+1]
    for (int i = 1; i < n - 1; ++i) {
        b[i] = do_ * x0[i - 1] + dp * x0[i] + do_ * x0[i + 1];
    }
    
    // Última fila: b[n-1] = do * x0[n-2] + dp2 * x0[n-1]
    b[n - 1] = do_ * x0[n - 2] + dp2 * x0[n - 1];
    
    return b;
}

// Realiza un paso ADI en la dirección x
std::vector<std::complex<float>> custom_solvers::adi_x(
        const std::vector<std::complex<float>>& phi,  // Campo complejo 2D aplanado (row-major)
        int Nx,                                       // Número de puntos en dirección x
        int Ny,                                       // Número de puntos en dirección y
        float eps,                                    // Valor pequeño para evitar división por cero
        float k,                                      // Número de onda
        float dz,                                     // Paso en dirección z
        float dx                                      // Paso en dirección x
) {
    // Crear vector para el resultado intermedio
    std::vector<std::complex<float>> phi_inter(phi.size());
    
    // Calcular factor ung (uso la notación compleja i en lugar de j)
    std::complex<float> ung(0.0f, dz / (4.0f * k * dx * dx));
    
    // Para cada fila en dirección y
    for (int j = 0; j < Ny; ++j) {
        // Extraer la columna j actual para procesarla
        std::vector<std::complex<float>> phi_col(Nx);
        for (int i = 0; i < Nx; ++i) {
            phi_col[i] = phi[i * Ny + j];  // Acceso row-major
        }
        
        // Calcular ratios para condiciones de contorno
        std::complex<float> ratio_x0, ratio_xn;
        
        // Ratio en el borde izquierdo
        if (std::abs(phi_col[1]) < eps) {
            ratio_x0 = std::complex<float>(1.0f, 0.0f);
        } else {
            ratio_x0 = phi_col[0] / phi_col[1];
        }
        
        // Ratio en el borde derecho
        if (std::abs(phi_col[Nx-2]) < eps) {
            ratio_xn = std::complex<float>(1.0f, 0.0f);
        } else {
            ratio_xn = phi_col[Nx-1] / phi_col[Nx-2];
        }
        
        // Calcular coeficientes para la matriz del lado derecho (B)
        std::complex<float> dp1_B = -2.0f * ung + std::complex<float>(1.0f, 0.0f) + ung * ratio_x0;
        std::complex<float> dp2_B = -2.0f * ung + std::complex<float>(1.0f, 0.0f) + ung * ratio_xn;
        std::complex<float> dp_B = -2.0f * ung + std::complex<float>(1.0f, 0.0f);
        std::complex<float> do_B = ung;
        
        // Calcular el vector b del lado derecho de la ecuación
        std::vector<std::complex<float>> b = compute_b_vector(dp_B, dp1_B, dp2_B, do_B, phi_col);
        
        // Calcular coeficientes para la matriz del lado izquierdo (A)
        std::complex<float> dp1_A = 2.0f * ung + std::complex<float>(1.0f, 0.0f) - ung * ratio_x0;
        std::complex<float> dp2_A = 2.0f * ung + std::complex<float>(1.0f, 0.0f) - ung * ratio_xn;
        std::complex<float> dp_A = 2.0f * ung + std::complex<float>(1.0f, 0.0f);
        std::complex<float> do_A = -ung;
        
        // Resolver el sistema tridiagonal
        std::vector<std::complex<float>> result_col(Nx);
        custom_thomas_solver(dp_A, dp1_A, dp2_A, do_A, b, result_col);
        
        // Guardar la solución en la matriz intermedia
        for (int i = 0; i < Nx; ++i) {
            phi_inter[i * Ny + j] = result_col[i];  // Acceso row-major
        }
    }
    
    return phi_inter;
}


// Realiza un paso ADI en la dirección y
std::vector<std::complex<float>> custom_solvers::adi_y(
        const std::vector<std::complex<float>>& phi,  // Campo complejo 2D aplanado (row-major)
        int Nx,                                       // Número de puntos en dirección x
        int Ny,                                       // Número de puntos en dirección y
        float eps,                                    // Valor pequeño para evitar división por cero
        float k,                                      // Número de onda
        float dz,                                     // Paso en dirección z
        float dy                                      // Paso en dirección y
) {
    // Crear vector para el resultado intermedio
    std::vector<std::complex<float>> phi_inter(phi.size());
    
    // Calcular factor ung (uso la notación compleja i en lugar de j)
    std::complex<float> ung(0.0f, dz / (4.0f * k * dy * dy));
    
    // Para cada fila en dirección x
    for (int i = 0; i < Nx; ++i) {
        // Extraer la fila i actual para procesarla
        std::vector<std::complex<float>> phi_row(Ny);
        for (int j = 0; j < Ny; ++j) {
            // Asumiendo que phi está almacenado en formato (x, y) -> índice = x*Ny + y
            phi_row[j] = phi[i*Ny + j];
        }
        
        // Calcular ratios para condiciones de contorno
        std::complex<float> ratio_y0, ratio_yn;
        
        // Ratio en el borde inferior
        if (std::abs(phi_row[1]) < eps) {
            ratio_y0 = std::complex<float>(1.0f, 0.0f);
        } else {
            ratio_y0 = phi_row[0] / phi_row[1];
        }
        
        // Ratio en el borde superior
        if (std::abs(phi_row[Ny-2]) < eps) {
            ratio_yn = std::complex<float>(1.0f, 0.0f);
        } else {
            ratio_yn = phi_row[Ny-1] / phi_row[Ny-2];
        }
        
        // Calcular coeficientes para la matriz del lado derecho (B)
        std::complex<float> dp1_B = -2.0f * ung + std::complex<float>(1.0f, 0.0f) + ung * ratio_y0;
        std::complex<float> dp2_B = -2.0f * ung + std::complex<float>(1.0f, 0.0f) + ung * ratio_yn;
        std::complex<float> dp_B = -2.0f * ung + std::complex<float>(1.0f, 0.0f);
        std::complex<float> do_B = ung;
        
        // Calcular el vector b del lado derecho de la ecuación
        std::vector<std::complex<float>> b = compute_b_vector(dp_B, dp1_B, dp2_B, do_B, phi_row);
        
        // Calcular coeficientes para la matriz del lado izquierdo (A)
        std::complex<float> dp1_A = 2.0f * ung + std::complex<float>(1.0f, 0.0f) - ung * ratio_y0;
        std::complex<float> dp2_A = 2.0f * ung + std::complex<float>(1.0f, 0.0f) - ung * ratio_yn;
        std::complex<float> dp_A = 2.0f * ung + std::complex<float>(1.0f, 0.0f);
        std::complex<float> do_A = -ung;
        
        // Resolver el sistema tridiagonal
        std::vector<std::complex<float>> result_row(Ny);
        custom_thomas_solver(dp_A, dp1_A, dp2_A, do_A, b, result_row);
        
        // Guardar la solución en la matriz intermedia
        for (int j = 0; j < Ny; ++j) {
            phi_inter[i*Ny + j] = result_row[j];
        }
    }
    
    return phi_inter;
}


// Aplica el operador no lineal (medio paso)
std::vector<std::complex<float>> custom_solvers::half_nonlinear(
        const std::vector<std::complex<float>>& phi,  // Campo complejo
        float k_sample,                               // Número de onda en el medio
        float n2_sample,                              // Índice de refracción no lineal
        float dz                                      // Paso en dirección z
) {
    int size = phi.size();
    std::vector<std::complex<float>> result(size);
    
    // Calcular el factor de fase no lineal para cada punto
    for (int i = 0; i < size; ++i) {
        // Calcular |phi|²
        float intensity = std::norm(phi[i]);  // std::norm devuelve |z|²
        
        // Calcular exp(i * k * n2 * dz/2 * |phi|²)
        std::complex<float> phase(0.0f, k_sample * n2_sample * dz/2.0f * intensity);
        std::complex<float> phase_factor = std::exp(phase);
        
        // Aplicar el factor de fase al campo
        result[i] = phase_factor * phi[i];
    }
    
    return result;
}

// Aplica la absorción lineal (medio paso)
std::vector<std::complex<float>> custom_solvers::half_linear_absorption(
        const std::vector<std::complex<float>>& phi,  // Campo complejo
        float alpha,                                  // Coeficiente de absorción lineal
        float dz                                      // Paso en dirección z
) {
    int size = phi.size();
    std::vector<std::complex<float>> result(size);
    
    // Calcular el factor de atenuación exp(-alpha * dz/4)
    float attenuation = std::exp(-alpha * dz/4.0f);
    
    // Aplicar la atenuación a cada punto del campo
    for (int i = 0; i < size; ++i) {
        result[i] = attenuation * phi[i];
    }
    
    return result;
}

// Aplica la absorción de dos fotones (medio paso)
std::vector<std::complex<float>> custom_solvers::half_2photon_absorption(
        const std::vector<std::complex<float>>& phi,  // Campo complejo
        float beta,                                   // Coeficiente de absorción de dos fotones
        float dz                                      // Paso en dirección z
) {
    int size = phi.size();
    std::vector<std::complex<float>> result(size);
    
    // Aplicar la atenuación no lineal a cada punto del campo
    for (int i = 0; i < size; ++i) {
        // Calcular |phi|²
        float intensity = std::norm(phi[i]);  // std::norm devuelve |z|²
        
        // Calcular exp(-beta * dz/4 * |phi|²)
        float attenuation = std::exp(-beta * dz/4.0f * intensity);
        
        // Aplicar la atenuación al campo
        result[i] = attenuation * phi[i];
    }
    
    return result;
}