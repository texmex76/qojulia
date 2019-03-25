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
m = 10
ωr = (k^2) / (2*m)

Δc = -10 / ωr
η = 50 / ωr
U0 = -1 / ωr

# Bases
b_position = PositionBasis(xmin, xmax, Nsteps)
b_fock = FockBasis(N_cutoff)

# Operators for atom
p = momentum(b_position)

# Operators for cavity
a = destroy(b_fock) ⊗ one(b_position)
ad = dagger(a)
H_kin = one(b_fock) ⊗ p^2/2m

dimension = N_cutoff * Nsteps + Nsteps
ψ0 = Ket(b_fock ⊗ b_position, rand(dimension) + im*rand(dimension))
α0 = rand(1)[1] * im*rand(1)[1]
κ = 10

potential = x -> cos(k*x)^2
H_cos_sqr = (one(b_fock) ⊗ potentialoperator(b_position, potential))

pump = x -> cos(k*x)
H_cos = (one(b_fock) ⊗ potentialoperator(b_position, pump))

id = one(b_fock) ⊗ identityoperator(b_position)

function propagate(α, ψ)
    dt_ψ = -im*(H_kin + abs2(α)*H_cos_sqr + 2*real(α)*η*H_cos + norm(ψ)^2*id)*ψ
    ψ += dt_ψ
    normalize!(ψ)
    dt_α = -im*((-Δc + U0*expect(H_cos_sqr, ψ) - im*κ)*α + η*expect(H_cos, ψ))
    α += dt_α
    return α, ψ
end

alphas = []
psis = []
alpha, psi = propagate(α0, ψ0)
push!(alphas, alpha)
push!(psis, psi)

for i in [1:1:100;]
    alpha, psi = propagate(alphas[size(alphas)[1]], psis[size(psis)[1]])
    push!(alphas, alpha)
    push!(psis, psi)
end

# println(abs(expect(H_cos, psis[100])))

psis = @. ptrace(psis, 1)

density = []
for i in [1:1:101;]
    push!(density, diag(psis[i].data))
end

norms = []
for i in [1:1:101;]
    push!(norms, norm(density[i]))
end

# clf()
# points = [1:1:301;]
# plot(points, norms)
# gcf()

clf()
points = [1:1:Nsteps;]
plot(points, abs2.(density[101]))
gcf()
