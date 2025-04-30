**Split-Step Beam Propagation Method (BPM) with ADI Crank–Nicolson**

This README outlines the mathematical formulation and numerical implementation of the split-step BPM using an Alternating-Direction Implicit (ADI) Crank–Nicolson scheme for solving the paraxial wave equation in two transverse dimensions.

---

## 1. Governing Equation

We solve the scalar **paraxial wave equation** for the complex field \(\Psi(x,y,z)\):

\[
\frac{\partial \Psi}{\partial z} =
\underbrace{\frac{i}{2k} \nabla_\perp^2 \Psi}_{\text{diffraction}}
+ \underbrace{i k n_2 |\Psi|^2 \Psi}_{\text{nonlinear phase}} \,.
  \]

- \(k = 2\pi / \lambda\): vacuum wavenumber
- \(\nabla_\perp^2 = \partial^2/\partial x^2 + \partial^2/\partial y^2\)
- \(n_2\): nonlinear refractive index coefficient


## 2. Operator Splitting (Split‑Step)

We split propagation over a small step \(\Delta z\) into two operators:

1. **Nonlinear operator**: \(\mathcal{N}[\Psi] = i k n_2 |\Psi|^2 \Psi\)
2. **Linear (diffraction) operator**: \(\mathcal{L}[\Psi] = \tfrac{i}{2k} \nabla_\perp^2 \Psi\)

Using a symmetric (Strang) split:

\[
\Psi(z + \Delta z) \approx e^{\tfrac{\Delta z}{2} \mathcal{N}} \, e^{\Delta z \mathcal{L}} \, e^{\tfrac{\Delta z}{2} \mathcal{N}} \, \Psi(z).
\]

- The **nonlinear phase step** is applied in closed form:
  \[ \Psi \to \Psi \exp\bigl(i k n_2 |\Psi|^2 \tfrac{\Delta z}{2}\bigr). \]

- The **linear diffraction** remains; we discretize and march it using ADI Crank–Nicolson.


## 3. Crank–Nicolson Discretization for Diffraction

Starting from:
\[ \frac{\partial \Psi}{\partial z} = \tfrac{i}{2k} \nabla_\perp^2 \Psi, \]
Crank–Nicolson gives for one full \(\Delta z\) step:

\[
\frac{\Psi^{n+1} - \Psi^n}{\Delta z}
= \frac{i}{4k} \nabla_\perp^2 \bigl(\Psi^{n+1} + \Psi^n\bigr).
\]

Rearrange into a matrix equation:

\[
\Bigl( I - \tfrac{i \Delta z}{4k} \nabla_\perp^2 \Bigr) \Psi^{n+1}
= \Bigl( I + \tfrac{i \Delta z}{4k} \nabla_\perp^2 \Bigr) \Psi^n.
\]

In two dimensions, directly inverting the full \(xy\)-Laplacian is costly. Instead, we use **ADI**:

1. **Half-step in x**:
   \[
   \Bigl( I - \tfrac{i \Delta z}{4k} \partial_x^2 \Bigr) \Psi^{*}
   = \Bigl( I + \tfrac{i \Delta z}{4k} \partial_x^2 \Bigr) \Psi^n.
   \]
2. **Half-step in y**:
   \[
   \Bigl( I - \tfrac{i \Delta z}{4k} \partial_y^2 \Bigr) \Psi^{n+1}
   = \Bigl( I + \tfrac{i \Delta z}{4k} \partial_y^2 \Bigr) \Psi^{*}.
   \]

Each sub‑step involves solving a **tridiagonal** system along one coordinate, which is efficiently done via the Thomas algorithm.


## 4. Spatial Discretization

On a uniform grid:

- Grid points \(x_i = (i - N_x/2) \Delta x\), \(i=0,\dots,N_x-1\)
- Second derivative via central differences:
  \[
  \partial_x^2 \Psi_i \approx \frac{\Psi_{i+1} - 2\Psi_i + \Psi_{i-1}}{(\Delta x)^2}.
  \]

Thus each ADI half‑step leads to a banded (tridiagonal) linear system:

\[
\Bigl(1 + r\Bigr) \Psi_{i}^{*} - \tfrac{r}{2} \bigl(\Psi_{i+1}^{*} + \Psi_{i-1}^{*}\bigr)
= \Bigl(1 - r\Bigr) \Psi_{i}^{n} + \tfrac{r}{2} \bigl(\Psi_{i+1}^{n} + \Psi_{i-1}^{n}\bigr),
\]

with \(r = i\,\Delta z/(4k\Delta x^2)\). A similar system applies in \(y\).


## 5. Full Propagation Loop

1. **Initialize** \(\Psi(x,y,0)\) at \(z=0\).
2. **For** each \(n=0,1,\dots,N_z-1\):
    1. Nonlinear half-step:  
       \(\Psi \leftarrow \Psi \exp\bigl(i k n_2 |\Psi|^2 \tfrac{\Delta z}{2}\bigr)\)
    2. Linear ADI (x then y): solve tridiagonal systems
    3. Nonlinear half-step again
3. **Output** \(\Psi(x,y,z_{\mathrm{final}})\) or intensity \( |\Psi|^2 \).


## 6. References

- Agrawal, G. P., *Nonlinear Fiber Optics*, 5th ed.
- Hardin, R. H. & Tappert, F. D., "Applications of the split-step Fourier method to the numerical solution of nonlinear and variable coefficient wave equations."
- Muir, T., "ADI Schemes for Diffraction Problems."

---


