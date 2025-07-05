#include <cstdio>
#include <complex>
#include <vector>
#include <iostream>
#include "LaserBeam.h"
#include "PrintHelpers.h"
#include "custom_solvers.h"


void create_beam_test1();
void thomas1_test();

int main() {
    thomas1_test();
    return 0;
}


void create_beam_test1(){
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
}

void thomas1_test(){
    using namespace std;
    using complexf = complex<float>;

    int N = 4;
    vector<complexf> a = {0.0f, 1.0f, 1.0f, 1.0f}; // a[0] is unused
    vector<complexf> b = {4.0f, 4.0f, 4.0f, 4.0f};
    vector<complexf> c = {1.0f, 1.0f, 1.0f, 0.0f}; // c[N-1] is unused
    vector<complexf> d = {complexf(5,1), complexf(6,2), complexf(6,1), complexf(5,0)};
    vector<complexf> x(N);

    custom_solvers::custom_thomas_solver(a, b, c, d, x);

    cout << "Solution x:\n";
    for (const auto& val : x)
        cout << val << '\n';
}




