using QuantumOptics
using PyPlot
using Printf

N_cutoff = 1

xmin = -3
xmax = 3
Nsteps = 200
m = .5
k = 1.6
Δc = 10
Δa = -1
η = 10
g0 = 10

# Bases
b_position = PositionBasis(xmin, xmax, Nsteps)
b_fock = FockBasis(N_cutoff)

# Operators for atom
p = momentum(b_position)

# Operators for cavity
a = destroy(b_fock) ⊗ one(b_position)
ad = dagger(a)

pot_cos_sqr = x -> g0^2/Δa*cos(k*x)^2
pot_cos = x -> -η*cos(k*x)
H_int = (one(b_fock) ⊗ potentialoperator(b_position, pot_cos_sqr))*ad*a
H_kin = one(b_fock) ⊗ p^2/2m
H_pump = (one(b_fock) ⊗ potentialoperator(b_position, pot_cos)) * (a + ad)
H_atom = -Δc*ad*a
H = H_kin + dense(H_int) + dense(H_pump) + H_atom

E, states = eigenstates((H + dagger(H))/2, 1);

# States has dimension N_cutoff*Nsteps + Nsteps
# We need to extract the points which only belong to the position basis
function extract_pos(Nsteps, N_cutoff, state)
    pos = []
    for i in 1:Nsteps
        index = (N_cutoff + 1)i
        push!(pos, state.data[index])
    end
    return pos
end

spts = extract_pos(Nsteps, N_cutoff, states[1])

xpoints = samplepoints(b_position)
pot = []
for i in xpoints
    push!(pot, -η*cos(k*i) + g0^2/Δa*cos(k*i)^2)
end

fig = figure()
host = fig.add_subplot(111)

par1 = host.twinx()

host.set_xlim(xmin, xmax)
host.set_ylim(-0.0014, 0.0297)
par1.set_ylim(-115.5123, 5.7583)

title("Transversal Pump with Operators")
host.set_xlabel(L"x")
host.set_ylabel(L"|\psi(x)|^2")
par1.set_ylabel(L"g_0 \Omega/\Delta a \cos(ki) + g_0^2/\Delta a \cos(ki)^2")

color1 = PyPlot.cm.Blues(200)
color2 = PyPlot.cm.Oranges(0.5)

host.plot(xpoints, abs2.(spts), color=color1)
par1.plot(xpoints, pot, color=color2, linestyle="--")

host.yaxis.label.set_color(color1)
par1.yaxis.label.set_color(color2)

gcf()

# println("host.set_ylim(" *@sprintf("%.4f", host.get_ylim()[1]) *", " *@sprintf("%.4f", host.get_ylim()[2]) *")")
# println("par1.set_ylim(" *@sprintf("%.4f", par1.get_ylim()[1]) *", " *@sprintf("%.4f", par1.get_ylim()[2]) *")")
