using QuantumOptics
using PyPlot

# System Parameters
m = 1
U0 = -1
k = .15
ωr = k^2/(2m)

# Position basis
xmin = -5
xmax = 5
Npoints = 100
b_position = PositionBasis(xmin, xmax, Npoints)
b_momentum = MomentumBasis(b_position)

# Transforms a state multiplied from the right side from real space
# to momentum space
T_px = particle.transform(b_momentum, b_position)

T_xp = dagger(T_px)
# x = position(b_position)
p = momentum(b_momentum)
id = identityoperator(b_position)

H_kin = LazyProduct(T_xp, p^2/2m, T_px)
H_cav = U0*id # ∝ cos(kx)^2

function fquantum(t, psi, u) # psi is quantum, u the classical part
    x = u[1]
    return LazySum(H_kin, H_cav*cos(x)^2)
end

function fclassical(t, psi, u, du)
    du[1] = 2ωr*u[2]
    du[2] = -2*U0*k*sin(k*u[1])*cos(k*u[1])*expect(id, psi)
end

# Initial state
x0 = 1.5
p0 = 0
sigma0 = 0.6
ψ0 = gaussianstate(b_position, x0, p0, sigma0)
u0 = complex.([x0, p0])

ψsc0 = semiclassical.State(ψ0, u0)

T = [0:0.1:400;]
tout, ρt = semiclassical.schroedinger_dynamic(T, ψsc0, fquantum, fclassical)

pos = []
mom = []

for r=ρt
    push!(pos, real(r.classical[1]))
    push!(mom, real(r.classical[2]))
end

pot = []
for i in pos
    push!(pot, U0*cos(k*i)^2)
end

figure(figsize=(10, 5))
subplot(221)
plot(tout, pos)
xlabel(L"t")
ylabel(L"x")

subplot(222)
plot(tout, mom.^2)
xlabel(L"t")
ylabel(L"E_{kin}")

subplot(223)
plot(tout, pot)
xlabel(L"t")
ylabel(L"U_0 \cos(kx)^2")

tight_layout()
gcf()
