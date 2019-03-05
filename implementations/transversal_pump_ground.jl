using QuantumOptics
using PyPlot

xmin = -5
xmax = 5
Nsteps = 200
m = 0.5
U0 = 10
k = 1

b_position = PositionBasis(xmin, xmax, Nsteps)
xpoints = samplepoints(b_position)

x = position(b_position)
p = momentum(b_position)

potential = x -> U0*cos(k*x)^2 + U0*cos(k*x)
V = potentialoperator(b_position, potential)
Hkin = p^2/2m
H = Hkin + dense(V)

E, states = eigenstates((H + dagger(H))/2, 1);

pot = []
for i in xpoints
    push!(pot, U0*cos(k*i)^2 + U0*cos(k*i))
end

figure(figsize=(8, 3))

subplot(221)
plot(xpoints, abs2.(states[1].data))
xlabel(L"x")
ylabel(L"|\psi(x)|^2")

subplot(222)
plot(xpoints, pot, "C1")
xlabel(L"x")
ylabel(L"V")

gcf()
