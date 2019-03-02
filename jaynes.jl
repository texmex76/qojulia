using QuantumOptics
using PyPlot

# Parameters
N_cutoff = 10

ωc = 0.1
ωa = 0.1
Ω = 1.

# Bases
b_fock = FockBasis(N_cutoff)
b_spin = SpinBasis(1//2)
b = b_fock ⊗ b_spin

# Fundamental operators
a = destroy(b_fock)
at = create(b_fock)
n = number(b_fock)

sm = sigmam(b_spin)
sp = sigmap(b_spin)
sz = sigmaz(b_spin)

# Hamiltonian
Hatom = ωa*sz/2
Hfield = ωc*n
Hint = Ω*(at⊗sm + a⊗sp)
H = one(b_fock)⊗Hatom + Hfield⊗one(b_spin) + Hint

α = 1.
ψ0 = coherentstate(b_fock, α) ⊗ spindown(b_spin)

# Integration time
T = [0:0.1:20;]

# Schroedinger time evolution
tout, ψt = timeevolution.schroedinger(T, ψ0, H)

exp_n = real(expect(n ⊗ one(b_spin), ψt))
exp_sz = real(expect(one(b_fock) ⊗ sz, ψt))

figure(figsize=(9,3))
subplot(1,2,1)
ylim([0, 2])
plot(T, exp_n);
xlabel(L"T")
ylabel(L"\langle n \rangle")

subplot(1,2,2)
ylim([-1, 1])
plot(T, exp_sz);
xlabel(L"T")
ylabel(L"\langle \sigma_z \rangle")

tight_layout()
gcf()
