# Pan-S7-Code
Repository contains the original FORTRAN code routine to generate the Pan S7 Rayleigh-Brillouin scattering lineshape for both spontaneous and coherent Rayleigh-Brillouin scattering based on the linearized Wang-Chang–
Uhlenbeck equation in kinetic theory. Contains the following files:

- crbs_mol7_1.for - Main code.

- lubksb.for - Subroutine to solve the set of linear equations AX=B.

- ludcmp.for - Subroutine to replace a given matrix by the LU decomposition.

- qsimp.for - Subroutine to determine w0, the plasma dispersion function.

- trapzd.for - Called in qsimp for numerical integration.

Also contained is a Matlab version of the Pan S7 FORTRAN code originally written by Emmanuel Stockman when he was a graduate student at Princeton University for those not savvy to FORTRAN.

The FORTRAN code was originally written by Mikhail Shneider and Xingguo Pan, c. 2004.
The theory is described in the following reference:

X. Pan, M. N. Shneider, and R. B. Miles, “Coherent Rayleigh-Brillouin scattering in molecular gases,” Phys. Rev. A, vol. 69, no. 3, p. 033814, Mar. 2004.

Any time this code (or transcript of it) is used in work related to an academic paper or publication, it MUST be mentioned that the S7 code was originally developed by Princeton and referenced by the above citation.


David Feng, 2017/9/7
