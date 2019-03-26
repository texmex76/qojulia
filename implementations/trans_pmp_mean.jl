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

κ = .9*ωr
Δc = -1 / ωr
# η = 20 / ωr
U0 = -1 / ωr
gc = 10*ωr*λ
# gc = 1

function evaluate(η)
η = η / ωr

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
# psis = @. ptrace(psis, 1)
#
# density = []
# for i in [1:1:iterations + 1;]
#     push!(density, diag(psis[i].data))
# end
#
# norms = []
# for i in [1:1:iterations + 1;]
#     push!(norms, norm(density[i]))
# end
#
# # clf()
# # points = [1:1:301;]
# # plot(points, norms)
# # gcf()
#
# clf()
# points = [1:1:Nsteps;]
# plot(points, abs2.(density[iterations + 1]))
# gcf()
# #### end visualize density


θ = abs(expect(H_cos, psis[iterations + 1]))
return θ
end

points = [0:1:15;]
θ = @. evaluate(points)
figure()
clf()
suptitle("Transversal Pump Order Parameter")
xlabel(L"\eta")
ylabel(L"\Theta = \langle \cos(kx) \rangle")
plot(points, θ)
gcf()
