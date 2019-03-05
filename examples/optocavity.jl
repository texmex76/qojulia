using QuantumOptics
using PyPlot

# Parameters
ω_mech = 10.
Δ = -ω_mech

# Constants
g = 1.
η = 2.
κ = 1.

# Basis
b_cav = FockBasis(4)
b_mech = FockBasis(10)

# Operators Cavity
a = destroy(b_cav) ⊗ one(b_mech)
at = create(b_cav) ⊗ one(b_mech)

# Operators Oscillator
b = one(b_cav) ⊗ destroy(b_mech)
bt = one(b_cav) ⊗ create(b_mech)

# Hamilton Operator
H_cav = -Δ*at*a + η*(a+at)
H_mech = ω_mech*bt*b
H_int = -g*(bt+b)*at*a

H = H_cav + H_mech + H_int

J = [a]
rates = [κ]

ψ0 = fockstate(b_cav,0) ⊗ fockstate(b_mech,2)

T = [0:0.2:50;]
tout, ρt = timeevolution.master(T,ψ0,H,J;rates=rates)

figure(figsize=(9, 3))

subplot(1, 2, 1)
title("Mechanical Oscillator")
plot(T, real(expect(bt*b, ρt)))
xlabel(L"t \kappa")
ylabel(L"\langle n_m \rangle")

subplot(1, 2, 2)
title("Cavity")
plot(T, real(expect(at*a, ρt)))
xlabel(L"t \kappa")
ylabel(L"\langle n_c \rangle")

tight_layout()
gcf()
