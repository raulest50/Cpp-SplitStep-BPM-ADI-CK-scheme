# Split-Step Beam Propagation Method (BPM) with ADI Crank–Nicolson

This README outlines the mathematical formulation and numerical implementation of the split-step BPM using an Alternating-Direction Implicit (ADI) Crank–Nicolson scheme for solving the paraxial wave equation in two transverse dimensions.

---

## 1. Governing Equation

We solve the scalar paraxial wave equation for the complex field $\Psi(x,y,z)$:

```math
\frac{\partial \Psi}{\partial z} = \frac{i}{2k} \nabla_\perp^2 \Psi + i k n_2 |\Psi|^2 \Psi
```

where:
- $k = 2\pi / \lambda$ (vacuum wavenumber)
- $\nabla_\perp^2 = \partial^2/\partial x^2 + \partial^2/\partial y^2$ (transverse Laplacian)
- $n_2$ = nonlinear refractive index coefficient

---

## 2. Operator Splitting (Split-Step)

Using Strang splitting over a small step $\Delta z$:

```math
\Psi(z + \Delta z) \approx e^{\frac{\Delta z}{2}N} \; e^{\Delta z L} \; e^{\frac{\Delta z}{2}N} \; \Psi(z)
```

- **Nonlinear operator**: $N[\Psi] = i k n_2 |\Psi|^2 \Psi$
- **Linear (diffraction) operator**: $L[\Psi] = \tfrac{i}{2k} \nabla_\perp^2 \Psi$

The nonlinear half-step is applied in closed form:

```math
\Psi \to \Psi \exp\Bigl(i k n_2 |\Psi|^2 \tfrac{\Delta z}{2}\Bigr)
```

---

## 3. Crank–Nicolson Discretization for Diffraction

Starting from:

```math
\frac{\partial \Psi}{\partial z} = \frac{i}{2k} \nabla_\perp^2 \Psi
```

Crank–Nicolson over a full $\Delta z$ gives:

```math
\frac{\Psi^{n+1} - \Psi^n}{\Delta z} = \frac{i}{4k} \nabla_\perp^2 \bigl(\Psi^{n+1} + \Psi^n\bigr)
```

Rearrange:

```math
\Bigl(I - \frac{i\,\Delta z}{4k} \nabla_\perp^2\Bigr)\Psi^{n+1} = \Bigl(I + \frac{i\,\Delta z}{4k} \nabla_\perp^2\Bigr)\Psi^n
```

Direct inversion of the 2D Laplacian is expensive. Instead, use **ADI**:

1. **Half-step in $x$**:

   $$\Bigl(I - \frac{i\,\Delta z}{4k}\,\partial_x^2\Bigr)\Psi^* = \Bigl(I + \frac{i\,\Delta z}{4k}\,\partial_x^2\Bigr)\Psi^n$$

2. **Half-step in $y$**:

   $$\Bigl(I - \frac{i\,\Delta z}{4k}\,\partial_y^2\Bigr)\Psi^{n+1} = \Bigl(I + \frac{i\,\Delta z}{4k}\,\partial_y^2\Bigr)\Psi^*$$

Each sub-step is a tridiagonal system solved via the Thomas algorithm.

---

## 4. Spatial Discretization

On a uniform grid:


$$x_i = \bigl(i - \tfrac{N_x}{2}\bigr)\Delta x,  \quad  y_j = \bigl(j - \tfrac{N_y}{2}\bigr)\Delta y$$


Central differences for second derivatives:

```math
\frac{\partial^2 \Psi}{\partial x^2}\Big|_i \approx \frac{\Psi_{i+1} - 2\Psi_i + \Psi_{i-1}}{(\Delta x)^2}
```

Resulting tridiagonal equation for the $x$ half-step (similar for $y$):

```math
(1 + r)\Psi_i^* - \frac{r}{2}(\Psi_{i+1}^* + \Psi_{i-1}^*) = (1 - r)\Psi_i^n + \frac{r}{2}(\Psi_{i+1}^n + \Psi_{i-1}^n)
```

with

```math
r = \frac{i\,\Delta z}{4\,k\,(\Delta x)^2}
```

---

## 5. Full Propagation Loop

1. **Initialize** $\Psi(x,y,0)$ on the $(x,y)$ grid.
2. **For** each step $n = 0,1,\dots,N_z-1$:
    1. Nonlinear half-step:
       
       $$\Psi \leftarrow \Psi \exp\Bigl(i k n_2 |\Psi|^2 \tfrac{\Delta z}{2}\Bigr)$$
       
    3. Linear ADI:
        - solve $x$–direction tridiagonal → $\Psi^*$
        - solve $y$–direction tridiagonal → $\Psi^{n+1}$
    4. Nonlinear half-step again
3. **Output** $\Psi(x,y,z_{\mathrm{final}})$ or intensity $|\Psi|^2$.

---

## 6. References

- Agrawal, G. P., *Nonlinear Fiber Optics*, 5th ed.
- Hardin & Tappert, “Split-step Fourier method for nonlinear wave equations.”
- Descriptions of the Thomas algorithm for tridiagonal systems.

