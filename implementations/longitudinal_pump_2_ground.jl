using QuantumOptics
using PyPlot

N_cutoff = 4

xmin = -3
xmax = 3
Nsteps = 200
m = 0.5
U0 = -1
k = 1
η = 1
Δc = -10

b_position = PositionBasis(xmin, xmax, Nsteps)
b_fock = FockBasis(N_cutoff)
xpoints = samplepoints(b_position)

# Atom
x = position(b_position)
p = momentum(b_position)

# Modes
a = destroy(b_fock) ⊗ one(b_position)
ad = dagger(a)

potential = x -> U0*cos(k*x)^2
V = one(b_fock) ⊗ potentialoperator(b_position, potential)
H_kin = one(b_fock) ⊗ p^2/2m
H_pump = η*(a + ad)
H_atom = Δc*ad*a
H = H_kin + dense(V)*ad*a + H_pump + H_atom

E, states = eigenstates((H + dagger(H))/2, 1);

# steps = N_cutoff * Nsteps + Nsteps
# print(steps)
# increment = (xmax - xmin)/steps
# xpts = [xmin:increment:xmax-increment;]
#
# pot = []
# for i in xpts
#     push!(pot, U0*cos(k*i)^2)
# end
#
# figure(figsize=(8, 3))
#
# subplot(221)
# plot(xpts, abs2.(states[1].data))
# xlabel(L"x")
# ylabel(L"|\psi(x)|^2")
#
# subplot(222)
# plot(xpts, pot, "C1")
# xlabel(L"x")
# ylabel(L"U_0\cos(kx)^2")
#
# gcf()

xpts = []
for i in 1:Nsteps
    push!(xpts, (N_cutoff + 1)i)
end

spts = []
for i in xpts
    push!(spts, states[1].data[i])
end

pot = []
for i in xpoints
    push!(pot, U0*cos(k*i)^2)
end

figure(figsize=(8, 3))

subplot(221)
plot(xpoints, abs2.(spts))
xlabel(L"x")
ylabel(L"|\psi(x)|^2")

subplot(222)
plot(xpoints, pot, "C1")
xlabel(L"x")
ylabel(L"U_0\cos(kx)^2")

gcf()
