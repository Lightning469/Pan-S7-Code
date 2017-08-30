# Pan-S7-Code
FORTRAN code to generate Pan S7 Rayleigh-Brillouin scattering lineshape. Contains the following files:

crbs_mol7_1.for - Main code.

lubksb.for - Subroutine to solve the set of linear equations AX=B.

ludcmp.for - Subroutine to replace a given matrix by the LU decomposition.

qsimp.for - Subroutine to determine w0, the plasma dispersion function.

trapzd.for - Called in qsimp for numerical integration.

Originally developed Mikhail Shneider and Xingguo Pan, c. 2004.
Theory is described in the following reference:

X. Pan, M. N. Shneider, and R. B. Miles, “Coherent Rayleigh-Brillouin scattering in molecular gases,” Phys. Rev. A, vol. 69, no. 3, p. 033814, Mar. 2004.
