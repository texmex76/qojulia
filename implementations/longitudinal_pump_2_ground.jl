using QuantumOptics
using PyPlot

N_cutoff = 1

xmin = -3
xmax = 3
Nsteps = 200
m = .5
U0 = 100
k = 1.3
η = 10
Δc = 1

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
H_kin = one(b_fock) ⊗ p^2/2m
H_pump = η*(a + ad)
H_atom = Δc*ad*a
H = H_kin + dense(H_int) + H_pump - H_atom

E, states = eigenstates((H + dagger(H))/2, 1);

# States has dimension N_cutoff*Nsteps + Nsteps
# We need to extract the points which only belong to the position basis
spts = []
for i in 1:Nsteps
    index = (N_cutoff + 1)i
    push!(spts, states[1].data[index])
end

xpoints = samplepoints(b_position)
pot = []
for i in xpoints
    push!(pot, U0*cos(k*i)^2)
end

fig = figure()
host = fig.add_subplot(111)

par1 = host.twinx()

host.set_xlim(xmin, xmax)

host.set_xlabel(L"x")
host.set_ylabel(L"|\psi(x)|^2")
par1.set_ylabel(L"U_0\cos(kx)^2")

color1 = PyPlot.cm.Blues(200)
color2 = PyPlot.cm.Oranges(0.5)

host.plot(xpoints, abs2.(spts), color=color1)
par1.plot(xpoints, pot, color=color2)

host.yaxis.label.set_color(color1)
par1.yaxis.label.set_color(color2)

gcf()
