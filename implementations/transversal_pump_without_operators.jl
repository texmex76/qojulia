using QuantumOptics
using PyPlot

xmin = -5
xmax = 5
Nsteps = 200
m = .5
k = 1.5
g0 = sqrt(10)
Δa = -1
Ω = 10

b_position = PositionBasis(xmin, xmax, Nsteps)
xpoints = samplepoints(b_position)

x = position(b_position)
p = momentum(b_position)

potential = x -> -g0^2/Δa*cos(k*x)^2 - g0/Δa*Ω*cos(k*x)
V = potentialoperator(b_position, potential)
Hkin = p^2/2m
H = Hkin + dense(V)

E, states = eigenstates((H + dagger(H))/2, 1);

pot = []
for i in xpoints
    push!(pot, -g0^2/Δa*cos(k*i)^2 - g0/Δa*Ω*cos(k*i))
end

fig = figure()
host = fig.add_subplot(111)

par1 = host.twinx()

host.set_xlim(xmin, xmax)
host.set_ylim(-24.7846, 44.7850)
par1.set_ylim(-0.0014, 0.0300)

title("Transversal pump without operators")
host.set_xlabel(L"x")
host.set_ylabel(L"-g_0/\Delta a \cos(ki)^2 - g_0/\Delta a \Omega \cos(ki)")
par1.set_ylabel(L"|\psi(x)|^2")

host.plot(xpoints, pot, "C0")
par1.plot(xpoints, abs2.(states[1].data), "C1", linestyle="--")

host.yaxis.label.set_color("C0")
par1.yaxis.label.set_color("C1")

gcf()
