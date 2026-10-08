# Shared CSD normalization correction

The native group-integrated flux basis is `b_i=sqrt(2i+1) P_i`, with zero-based
polynomial degree. The stored unnormalized tensor is
`W_jkw=(1/2) integral(P_j P_k P_w du)`.

The Legendre derivative identity is
`P_i'=sum_(k<i, i-k odd) (2k+1) P_k`.
Consequently the integrated energy derivative contains
`C_i C_j C_k^2 C_w [1-(-1)^(i-k)] S_w W_jkw`, where native indices are one-based
and `C_k^2=2k-1`. The previous implementation had only one factor `C_k`.
Degrees zero and one are unaffected; higher energy test rows were incorrect.

Both reference and optimized matrix assemblers are corrected in SN/GN 1D/2D/3D.
The GN energy derivative acts as the identity on angular coefficient blocks.
Its independent fixture uses two angular blocks and non-diagonal streaming
matrices, with manufactured source `sum_a N_a v phi/dx_a - v d(S phi)/dE`.
Thus it tests the angular-block layout rather than assuming an SN ordinate result.

The SN fixtures use separately coded fixed six-node Gauss–Legendre integration
of polynomials, not the native tensor, polynomial or source helpers. GN shares
this independent energy quadrature and adds a separate angular-block derivation.
Constant and affine stopping, energy orders2/3/4, both coupling layouts, two
nonunit energy widths and unequal spatial widths are exercised. Each fixture
checks trial coefficients, spatial outgoing traces and lower-energy trace;
optimized matrix solves are checked separately. Dividing the independently
integrated tensor by `C_k` reconstructs the old defect: higher-order comparisons
must fail, while degree-one controls remain unchanged.

Spatial polynomial degree is zero in these manufactured controls. They do not
qualify adaptive closures, higher spatial degrees, the complete optimized sweep,
positivity, physical dose, DPN angular quadrature or a diffusion operator.
No straggling term has been enabled by this correction. Current native SN/GN
still perform a downward energy sweep; conservative intergroup diffusion needs
an implicit coupled-energy solve and additional boundary conditions.
