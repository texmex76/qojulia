using QuantumOptics
using PyPlot

# Parameters
Nc = 16
γ = 1
g = γ/2
κ = 0.5γ
ωr = .15γ
Δc = -γ
Δa = -2γ
η = γ
tmax = 400
tsteps = 10*tmax
dt = tmax / tsteps
tlist = [0:dt:tmax;]

# Hilbert space
bc = FockBasis(Nc)
ba = SpinBasis(1//2)

# Operators
a = destroy(bc) ⊗ one(ba)
ad = dagger(a)
sm = one(bc) ⊗ sigmam(ba)
sp = dagger(sm)

# Hamiltonian
#H0 = -Δc*ad*a - Δa*sp*sm + η*(a + ad)
H0 = -Δc*ad*a - Δa*sp*sm + η*(a + ad)
Hx = g*(a*sp + ad*sm) # ∝ cos(x)

# Jump operators
J = [sqrt(2κ)*a, sqrt(2γ)*sm]
Jdagger = map(dagger, J)

function fquantum(t, psi, u) # psi is quantum, u the classical part
    x = u[1]
    return H0 + Hx*cos(x), J, Jdagger
end

adsm = ad*sm # Define to avoid doing a multiplicaiton at every timestep
# du is a vector containing the increments of the classical variables
function fclassical(t, psi, u, du)
    du[1] = 2ωr*u[2]
    du[2] = 2g*sin(u[1])*real(expect(adsm, psi))
    #du[2] = g*sin(u[1])*expect((a*sp + ad*sm), psi)
end

x0 = sqrt(2) # Some arbitraty initial position
p0 = 7 # Some arbitraty initial momentum
u0 = complex.([x0, p0])

ψ0 = fockstate(bc, 0) ⊗ spindown(ba)

ψsc0 = semiclassical.State(ψ0, u0)

tout, ρt = semiclassical.master_dynamic(tlist, ψsc0, fquantum, fclassical)

x = []
p = []

for r=ρt
    push!(x, real(r.classical[1]))
    push!(p, real(r.classical[2]))
end

n = real(expect(ad*a, ρt))
pe = real(expect(sp*sm, ρt))

figure(figsize=(10, 5))
subplot(221)
plot(tout, x)
xlabel(L"t")
ylabel(L"x")

subplot(222)
plot(tout, p.^2)
xlabel(L"t")
ylabel(L"E_{kin}")

subplot(223)
plot(tout, n)
xlabel(L"t")
ylabel(L"n")

subplot(224)
plot(tout, pe)
xlabel(L"t")
ylabel(L"\langle \sigma^+\sigma^- \rangle")

tight_layout()
gcf()
savefig("rk.svg")
