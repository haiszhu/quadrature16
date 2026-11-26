# quadrature16

tools for 16th order quadrature

## Layout
- `src/`: Fortran quadrature16 implementation and MATLAB MEX gateway.

- `utils/`: MATLAB quadrature helpers, which includes kernel and geometry routines borrowed from BIE3D and Helsing's myGMRES Ax=b solve.

- `examples/`: convergence test pipeline for part of the kernels.

## Supported kernels

- A set of kernels: $\mathcal{K}=\lbrace\Delta u=0;\ (\Delta+k^2)u=0;\ \mu\Delta\mathbf{u}-\nabla p=0,\ \nabla \cdot\mathbf{u}=0\rbrace\times\lbrace\mathcal{S},\mathcal{D},\mathcal{S}',\mathcal{D}'\rbrace\times\lbrace\text{on-surface},\text{off-surface}\rbrace$

- quadrature16: $\mathcal{K}_{\mathrm{quadrature16}}=\mathcal{K}$, no Stokes pressure

- wish list: complexification and better close evaluation

## Purpose

A cleaned-up collection of quadrature routines I wrote for another project, where the project needs careful handling of the underlying quadrature. $order=16$ only.

## To do list

* (to do) 