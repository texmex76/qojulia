using QuantumOptics
using PyPlot

N_cutoff = 1

xmin = -.5
xmax = .5
Nsteps = 16
m = .5
g0 = 10
Δa = -1
k = 2 * π
Δc = -10
η = 100
U0 = -1

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
pump = x -> η*cos(k*x)
H_pump = (one(b_fock) ⊗ potentialoperator(b_position, pump)) * (a + ad)
H_atom = -Δc*ad*a
H = H_kin + dense(H_int) + H_pump + H_atom

E, states = eigenstates((H + dagger(H))/2, 2);
println(E)

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

combined_1 = Ket[(states[1] + states[2]) / sqrt(2)]
combined_2 = Ket[(states[1] - states[2]) / sqrt(2)]

comb_states_1 = extract_pos(Nsteps, N_cutoff, combined_1[1])
comb_states_2 = extract_pos(Nsteps, N_cutoff, combined_2[1])

xpoints = samplepoints(b_position)
pot = []
for i in xpoints
    push!(pot, η*cos(k*i) + g0^2/Δa*cos(k*i)^2)
end

# fig = figure()
# host = fig.add_subplot(111)
#
# par1 = host.twinx()
# par2 = host.twinx()
#
# host.set_xlim(xmin, xmax)
# host.set_ylim(-999.9979, 1199.9555)
# par1.set_ylim(-0.0011, 0.0231)
# par2.set_ylim(-0.0011, 0.0231)
#
# title("Transversal Pump with Operators")
# host.set_xlabel(L"x")
# host.set_ylabel(L"U_0\cos(kx)^2")
# par1.set_ylabel(L"|\psi(x)|^2")
# par2.set_ylabel(L"|\psi(x)|^2")
#
# # color1 = PyPlot.cm.Blues(200)
# # color2 = PyPlot.cm.Oranges(0.5)
#
# host.plot(xpoints, pot, "C0")
# par1.plot(xpoints, abs2.(comb_states_1), "C1", linestyle="--")
# par2.plot(xpoints, abs2.(comb_states_2), "C2", linestyle="-.")
#
# host.yaxis.label.set_color("C0")
# par1.yaxis.label.set_color("C1")
# par2.yaxis.label.set_color("C2")
#
# gcf()

# println("host.set_ylim(" *@sprintf("%.4f", host.get_ylim()[1]) *", " *@sprintf("%.4f", host.get_ylim()[2]) *")")
# println("par1.set_ylim(" *@sprintf("%.4f", par1.get_ylim()[1]) *", " *@sprintf("%.4f", par1.get_ylim()[2]) *")")
clf()
pos_dense = ptrace(states[1], 1)
density = diag(pos_dense.data)
plot(xpoints, abs2.(density))
# ylim(0, 1)
gcf()
