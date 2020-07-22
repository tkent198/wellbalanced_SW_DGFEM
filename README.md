# wellbalanced_SW_DGFEM

## Ensuring 'well-balanced' shallow water flows via a discontinuous Galerkin finite element method: issues at lowest order

The discontinuous Galerkin finite element method (DGFEM) developed by Rhebergen et al. (2008) offers a robust method for solving systems of nonconservative hyperbolic partial differential equations but, as we show here, does not satisfactorily deal with topography in shallow water flows at lowest order (so-called DG0, or equivalently finite volume). In particular, numerical solutions of the space-DG0 discretised one-dimensional shallow water equations over varying topography are not truly `well-balanced'. A numerical scheme is well-balanced if trivial steady states are satisfied in the numerical solution; in the case of the shallow water equations, initialised rest flow should remain at rest for all times. Whilst the free-surface height and momentum remain constant and zero, respectively, suggesting that the scheme is indeed well-balanced, the fluid depth and topography evolve in time. This is both undesirable and unphysical, leading to incorrect numerical solutions for the fluid depth, and is thus a concern from a predictive modelling perspective. We expose this unsatisfactory issue, both analytically and numerically, and indicate a solution that combines the DGFEM formulation for nonconservative products with a fast and stable well-balanced finite-volume method. This combined scheme bypasses the offending issue and successfully integrates nonconservative hyperbolic shallow water-type models with varying topography at lowest order. We briefly discuss implications for the definition of a well-balanced scheme, and highlight applications when higher-order schemes may not be desired, which give further value to our finding beyond its exposure alone.

## Getting started
### Add Language, versions, etc.
* MATLAB '9.4.0.813654 (R2018a)'

## Files overview
### MATLAB

File name                   |  Summary
:--------------------------:|:--------------------------:
```...```          |  ...

---
![live_panel](MATLAB/mov/wetro2_Nk_100_Tend_1000.gif)
---