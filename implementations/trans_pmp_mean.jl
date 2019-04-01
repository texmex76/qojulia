using QuantumOptics
using PyPlot
using Printf
using LinearAlgebra
pygui(false)

N_c = 16
xmin = 0
xmax = 1
Nsteps = 32

k = 2*π
m = 10
ωr = (k^2) / (2*m)
λ = 2*π / k

κ = ωr
Δc = -1 / ωr
η = 30 / ωr
U0 = -1 / ωr
gc = 10*ωr*λ

# function evaluate(η)
# η = η / ωr

# Bases
b_position = PositionBasis(xmin, xmax, Nsteps)
b_fock = FockBasis(N_c)

# Operators for atom
p = momentum(b_position)

# Operators for cavity
a = destroy(b_fock) ⊗ one(b_position)
ad = dagger(a)
H_kin = (one(b_fock) ⊗ p^2/2m) / k^2

dimension = N_c * Nsteps + Nsteps
ψ0 = Ket(b_fock ⊗ b_position, rand(dimension) + im*rand(dimension))
α0 = rand(1)[1] * im*rand(1)[1]

potential = x -> cos(k*x)^2
H_cos_sqr = (one(b_fock) ⊗ potentialoperator(b_position, potential))

pump = x -> cos(k*x)
H_cos = (one(b_fock) ⊗ potentialoperator(b_position, pump))

id = one(b_fock) ⊗ identityoperator(b_position)

function propagate(α, ψ)
    dt_ψ = -im*(-H_kin + U0*abs2(α)*id + 2*real(α)*η*H_cos + norm(ψ)^2*id)*ψ
    ψ += dt_ψ
    normalize!(ψ)
    dt_α = -im*((-Δc + U0*expect(H_cos_sqr, ψ) - im*κ)*α + im*η + η*expect(H_cos, ψ))
    α += dt_α
    return α, ψ
end

alphas = []
psis = []
alpha, psi = propagate(α0, ψ0)
push!(alphas, alpha)
push!(psis, psi)

iterations = 200
for i in [1:1:iterations;]
    alpha, psi = propagate(alphas[size(alphas)[1]], psis[size(psis)[1]])
    push!(alphas, alpha)
    push!(psis, psi)
end

# #### begin visualize density
psis = @. ptrace(psis, 1)

density = []
for i in [1:1:iterations + 1;]
    push!(density, diag(psis[i].data))
end

xpoints = samplepoints(b_position)
pot = @. U0*cos(k*xpoints)^2*abs2(alphas[iterations + 1]) + η*cos(k*xpoints)*2*real(alphas[iterations + 1])

fig = figure()
host = fig.add_subplot(111)
par1 = host.twinx()

host.set_xlim(xmin, xmax)
host.set_ylim(-0.0013, 0.0278)
par1.set_ylim(-883.033083, 548.013803)

title("Transversal Pump Mean Field")
host.set_xlabel(L"x")
host.set_ylabel(L"|\psi(x)|^2")
par1.set_ylabel(L"U_0\cos(kx)^2\langle a^\dagger a\rangle + \eta \cos(kx) \langle a + a^\dagger \rangle")


host.plot(xpoints, abs2.(density[iterations + 1] ./ sqrt(Nsteps)), "C0")
par1.plot(xpoints, pot, "C1", linestyle="--")

host.yaxis.label.set_color("C0")
par1.yaxis.label.set_color("C1")

gcf()

# println("host.set_ylim(" *@sprintf("%.4f", host.get_ylim()[1]) *", " *@sprintf("%.4f", host.get_ylim()[2]) *")")
# println("par1.set_ylim(" *@sprintf("%.6f", par1.get_ylim()[1]) *", " *@sprintf("%.6f", par1.get_ylim()[2]) *")")

### end visualize density

# ## begin calculate order parameter
# θ = abs(expect(H_cos, psis[iterations + 1]))
# return θ
# end
#
# points = [0:1:15;]
# θ = @. evaluate(points)
# figure()
# clf()
# suptitle("Transversal Pump Order Parameter")
# xlabel(L"\eta")
# ylabel(L"\Theta = \langle \cos(kx) \rangle")
# plot(points, θ)
# gcf()
# ## end calculate order parameter
