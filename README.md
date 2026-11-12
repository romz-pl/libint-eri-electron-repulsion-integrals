# libint: Electron Repulsion Integrals (ERI) library

**Libint** is a C++ library for computing molecular integrals over Gaussian basis functions. Its best-known capability is the **electron repulsion integral (ERI)**, the most expensive part of many quantum chemistry calculations. It was developed mainly by Edward F. Valeev and collaborators, and it is widely used as a computational kernel inside larger electronic structure codes rather than as a standalone program.

## The integrals it computes

In Hartree–Fock, density functional theory, and post-Hartree–Fock methods, the basic objects are integrals over contracted Gaussian-type orbitals (GTOs). For a set of basis functions $\chi_\mu$, the two-electron integrals are

$$
(\mu\nu|\lambda\sigma) = \int\!\!\int \chi_\mu^*(\mathbf{r}_1)\,\chi_\nu(\mathbf{r}_1)\,\frac{1}{r_{12}}\,\chi_\lambda^*(\mathbf{r}_2)\,\chi_\sigma(\mathbf{r}_2)\, d\mathbf{r}_1\, d\mathbf{r}_2
$$

Libint evaluates these for contracted Cartesian Gaussians of arbitrary angular momentum, grouped into shells. It also handles the other one-electron and two-electron operators that electronic structure codes need:

- Overlap, kinetic energy, and nuclear attraction integrals
- Multipole moment integrals
- Electron repulsion integrals and their derivatives, which are needed for analytic gradients and Hessians
- Integrals for explicitly correlated methods, such as Gaussian geminal integrals used in F12 theory

Derivative integrals are particularly important because analytic gradients make geometry optimization and molecular dynamics far cheaper than finite differencing.

## Algorithms

Libint's core approach is the Obara–Saika scheme with Head-Gordon–Pople (HGP) recurrence relations. These reduce integrals over high angular momentum functions to a small set of auxiliary integrals, using horizontal and vertical recurrences that avoid the redundant work of naive expansions. The auxiliary quantities depend on the Boys function $F_m(T)$, which Libint evaluates with precomputed tables and asymptotic expansions. Rys quadrature is another standard route to the same integrals, and related codes use it extensively, but Libint's principal method is the recurrence-based one.

## Code generation

The distinguishing feature of Libint is its code generator. Integral classes are defined by the angular momenta of the four shells in a quartet, such as $(ps|dd)$. Instead of writing a hand-tuned kernel for each class, Libint uses a generator that emits specialized C++ for each class, with the recurrence steps expanded and optimized. The generated code is unrolled and organized to help the compiler vectorize it with SIMD instructions.

This approach matters because the number of possible shell combinations grows quickly with angular momentum. The generator makes it practical to support high angular momentum, with the maximum set at build time, and to add new operators or derivative orders without rewriting the kernels by hand.

## The C++ interface

Libint 2 uses C++ templates and an object-oriented interface. The main pieces are:

- **Shell** objects, which hold a center, the angular momentum, the primitive exponents, and the contraction coefficients.
- **BasisSet** objects, which collect shells for a molecule and a basis set such as cc-pVDZ.
- **Engine** objects, which are configured for a particular operator (Coulomb, overlap, kinetic, nuclear attraction, and so on), a maximum number of primitives, and a maximum angular momentum. Derivative order is also selected here.

A typical workflow is to create an Engine, loop over shell quartets, call `compute()`, and read the result from a buffer that Libint owns and reuses. Engines are not shared between threads, so parallel codes give each thread its own Engine. This design lets an application control memory, ordering, and parallel decomposition.

## Screening and performance

Most of the cost in an ERI calculation is avoidable, because many shell quartets contribute negligibly. Libint supports the standard Schwarz screening approach, which bounds each quartet using the diagonal integrals $(\mu\nu|\mu\nu)$. Codes that use Libint typically combine this with density-based screening so that quartets that cannot affect the final energy are skipped. Libint's kernels also support primitive-level screening and contraction-aware evaluation, which further reduce work.

## Where it is used

Libint is a dependency of several quantum chemistry packages. Psi4 is among the best known. The MPQC (Massively Parallel Quantum Chemistry) program, which Valeev's group developed, is built around it. Other codes use different ERI libraries. For example, PySCF relies on libcint, and several older or commercial programs use their own Rys-quadrature or Obara–Saika implementations.

## Licensing and practical notes

Libint is open source under the LGPL, so it can be linked into a larger application under the terms of that license. Building it requires choosing the maximum angular momentum and the operators you need, since the generated code depends on those choices. Compile times and library size grow with the angular momentum limit, so builds are often tuned to the basis sets in use.

