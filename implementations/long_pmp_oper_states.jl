using QuantumOptics
using PyPlot
using Printf
using LinearAlgebra
pygui(false)

N_cutoff = 16
xmin = 0
xmax = 1
Nsteps = 32

k = 2*π
ωr = 1

Δc = -30 * ωr
η = 50 * ωr
U0 = -10 * ωr

# function evaluate(η)
# η = η / ωr
# Bases
b_position = PositionBasis(xmin, xmax, Nsteps)
b_fock = FockBasis(N_cutoff)

# Operators for atom
p = momentum(b_position)

# Operators for cavity
a = destroy(b_fock) ⊗ one(b_position)
ad = dagger(a)

potential = x -> U0*cos(k*x)^2
H_int = (one(b_fock) ⊗ potentialoperator(b_position, potential))*ad*a
H_kin = (one(b_fock) ⊗ p^2) / k^2
H_pump = η*(a + ad)
H_cavity = -Δc*ad*a
H = H_kin + dense(H_int) + H_pump + H_cavity

E, states = eigenstates((H + dagger(H))/2, 1)

pos_dense = ptrace(states[1], 1)
density = diag(pos_dense.data)

ada_exp = expect(ad*a, states[1])

xpoints = samplepoints(b_position)
# pot = @. U0*cos(k*xpoints)^2*real(ada_exp)
#
# fig = figure()
# host = fig.add_subplot(111)
# par1 = host.twinx()
#
# host.set_xlim(xmin, xmax)
# host.set_ylim(-0.0001, 0.0017)
# par1.set_ylim(-299.233685, 14.249223)
#
# title(L"\eta = 50 \omega_r")
# host.set_xlabel(L"x")
# host.set_ylabel(L"|\psi(x)|^2")
# par1.set_ylabel(L"U_0\cos(kx)^2\langle a^\dagger a\rangle")
#
#
# host.plot(xpoints, abs2.(density ./ sqrt(Nsteps)), "C0")
# par1.plot(xpoints, pot, "C1", linestyle="--")
#
# host.yaxis.label.set_color("C0")
# par1.yaxis.label.set_color("C1")
# tight_layout()
#
# gcf()

# println("host.set_ylim(" *@sprintf("%.4f", host.get_ylim()[1]) *", " *@sprintf("%.4f", host.get_ylim()[2]) *")")
# println("par1.set_ylim(" *@sprintf("%.6f", par1.get_ylim()[1]) *", " *@sprintf("%.6f", par1.get_ylim()[2]) *")")

# pot_exp = x -> cos(k*x)^2
# cos_pot = one(b_fock) ⊗ potentialoperator(b_position, pot_exp)
# θ = real(expect(cos_pot, states[1]))
# return θ
# end
#
# points = [1:1:60;]
# θ = @. evaluate(points)
# figure()
# clf()
# suptitle("Longitudinal Pump Bunching Parameter")
# xlabel(L"\eta")
# ylabel(L"\langle \cos(kx)^2 \rangle")
# plot(points, θ)

b_momentum = MomentumBasis(b_position)
Tpx = transform(b_momentum, b_position)

states_p = Tpx * pos_dense
density_p = diag(states_p.data)

for item in density_p
    @printf("%.2f\n", abs(item))
end

clf()
bar([-16:1:15;], density_p, width=1)
gcf()
