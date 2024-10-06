
# Weakly first order melting of the 1/3 plateau in the Shastry-Sutherland model with iPEPS 

preprint: https://arxiv.org/pdf/2406.10689

The code simulates the thermal properties of the Shastry-Sutherland model (SSM) using iPEPS. The implementation of the U(1) symmetry is done using the iTensor library. The Shastry-Sutherland lattice is mapped onto the square lattice by considering the dimer basis, which allow us to compute the infinite contraction using the CTMRG algorithm.

The system is initialised at infinite temperature and then cooled down using the simple update scheme up to $\beta= 1/T$ using a $2 \times 2$ unit cell. The evaluation of local observables, and in particular the free energy can be done using a $2 \times 2$ or $6 \times 6$ unit cell. The infinite lattice is contracted using the $6 \times 6$ unit cell by projecting the boundary condition on the spin order of one of the ordered state. This in turn, improves the convergence of the CTMRG algorithm significantly in the ordered phase. In contrast, the infinite lattice is contracted using the 2x2 unit cell by using open boundary conditions.

We observe both the Z2 and Z3 symmetries of the ground state to restore at a unique temperature where the free energy has a kink, and conclude for the transition to occur in a single-step first-order transition.
