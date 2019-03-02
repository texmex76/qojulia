using QuantumOptics
using PyPlot

# System Parameters
m = 1
η = 1
k = 1
U0 = 1

N_cutoff = 10

# Position Basis
x_min = -5
x_max = 5
Npoints = 100
b_position = PositionBasis(x_min, x_max, Npoints)
b_momentum = MomentumBasis(b_position)

# Transforms a state multiplied from the right side from real space
# to momentum space
T_px = particle.transform(b_momentum, b_position)

T_xp = dagger(T_px)
x = position(b_position)
p = momentum(b_momentum)

# Fock basis
b_fock = FockBasis(N_cutoff)
a = destroy(b_fock)
at = create(b_fock)
n = number(b_fock)

H_kin = LazyProduct(T_xp, p^2/2m, T_px)
H_cav = η*cos(k*x)*(a+at)
H_int = U0*cos(k*x)^2*n

H = H_kin + H_cav + H_int
