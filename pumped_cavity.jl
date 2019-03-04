using QuantumOptics
using PyPlot

# Parameters
# Nc = 16
# κ = 20
# ωr = .15
# Δc = -1
# η = 1
# k = 1
# U0 = -1/8

Nc = 1
k = 1
m = 1
κ = 1
ωr = k^2/(2m)
κ = 200*ωr
Δc = -300*ωr
η = 50*ωr
U0 = -100*ωr

# Hilbert space
bc = FockBasis(Nc)

# Operators
a = destroy(bc)
ad = dagger(a)

# Hamiltonian
#H0 = -Δc*ad*a
H_cav = η*(a + ad) # ∝ cos(x)
H_int = U0*ad*a # ∝ cos(x)^2

# Jump operators
J = [sqrt(2κ)*a]
Jdagger = map(dagger, J)

function fquantum(t, psi, u) # psi is quantum, u the classical part
    x = u[1]
    return H_cav*cos(k*x) + H_int*cos(k*x)^2, J, Jdagger
end

ada = ad*a # Define to avoid doing a multiplicaiton at every timestep
# du is a vector containing the increments of the classical variables
function fclassical(t, psi, u, du)
    du[1] = 2ωr*u[2]
    du[2] = -η*k*sin(u[1])*expect(a + ad, psi)-2*U0*k*sin(k*u[1])*cos(k*u[1])*expect(ada, psi)
    #du[2] = g*sin(u[1])*expect((a*sp + ad*sm), psi)
end

x0 = sqrt(2) # Some arbitraty initial position
p0 = 7 # Some arbitraty initial momentum
u0 = complex.([x0, p0])

ψ0 = fockstate(bc, 0)

ψsc0 = semiclassical.State(ψ0, u0)

T = [0:0.1:400;]
tout, ρt = semiclassical.master_dynamic(T, ψsc0, fquantum, fclassical)

x = []
p = []

for r=ρt
    push!(x, real(r.classical[1]))
    push!(p, real(r.classical[2]))
end

n = real(expect(ada, ρt))

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

# subplot(224)
# plot(tout, pe)
# xlabel(L"t")
# ylabel(L"\langle \sigma^+\sigma^- \rangle")

tight_layout()
gcf()
#savefig("rk.svg")
