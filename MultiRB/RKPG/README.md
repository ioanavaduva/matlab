### Rational Krylov Solver

This folder contains a number of versions of the Rational Krylov solver created. Below, each is described.

### Prerequisites

Require access to rktoolbox, which can be obtained from here: http://guettel.com/rktoolbox/download.html


### Running the tests

To solve the Lyapunov equation with rank 1 rhs, 2-sided projection, need:
- driverRKPG.m: set up the problem by changing the size (n), rhs, and poles, as well as tolerance tol and maximum number of iterations maxit.
- RKPG.m: this is the solver; can choose the basis from:
    - get_rk_basis.m: rational Krylov basis (only works with rank 1 rhs).
    - get_pk_basis.m: standard (polynomial) Krylov basis (only tested for rank 1 rhs).
    - get_ek_basis.m: extended Krylov basis (only tested for rank 1 rhs).

To solve the Lyapunov equation with rank 1 rhs, 1-sided projection (on the left), need:
- oneside_RKPGdriver.m: set up the problem by changing the size (n), rhs, and poles, as well as tolerance tol and maximum number of iterations maxit.
- oneside_RKPG.m: this is the solver; uses:
    - get_rk_basis.m: rational Krylov basis (only works with rank 1 rhs).

To solve the Lyapunov equation with rank >1 rhs, 2-sided projection, need:
- driverBlockRKPG.m: set up the problem by changing the size (n), rhs, and poles, as well as tolerance tol and maximum number of iterations maxit; in driverBlockRKPG.m can change how iterations are counted by changing the solver:
    - RKPGblock.m: this solver counts iterations corresponding to each new column added to the basis; 
    - RKPGblock2.m: this solver counts iterations corresponding to each new use of a shift;
    - both solvers give the same solution and use:
        - get_start_basis.m generates the starting basis of same dimension as rhs.
        - get_rk_block_basis.m generates the block rational Krylov basis by adding a new column to existing ones.

To compare the bases with rat_krylov.m from rktoolbox use:
- CompareBases.m for rank 1 rhs; can choose size (n), rhs (f), number of poles (m - for random, k - for roots denominator).
- CompareBlockBases.m for higher rank rhs; can choose size (n), rhs (f), number of poles (m).

To compute different bounds: 
- DKSbound.m: computes the error bound at each iteration for the rational Krylov basis; from Theorem 4.8 in DKS paper: https://pdfs.semanticscholar.org/83f9/bf47d1d6372727a2ee5fef4bd3129e1c8e4a.pdf; need to input problem size (n) and number of iterations required for RKPG to converge.
- extended_bound.m: computes the error bound at each iteration for the extended Kryov basis; from eq (4.1) in Druskin, Simoncini paper: https://link.springer.com/content/pdf/10.1007/s00211-011-0366-3.pdf; need to input problem size (n) and number of iterations required for RKPG to converge.
- polynomial_bound.m: computes the error bound at each iteration for the standard (polynomial) Kryov basis; from Prop 3.1 in Simoncini, Druskin paper: https://www.jstor.org/stable/25663151?seq=1#metadata_info_tab_contents.
- u_out_product.m to be used in RKPG.m to compute the rational Krylov residual bound from Corollary 2.5, Beckermann paper: http://math.univ-lille1.fr/~bbecker/abstract/rational_Galerkin_rev.pdf

All other codes and subfolders in this folder are tests and plots run by me and can be ignored.