#
# module for dynamic spin bodies
#

module SpinBodies

using Statistics

export SpinBody, SpinLattice
export visualize, update_energies!, update_energies_fast!, update_lattice!, metropolis_step!
export prb, amp

mutable struct SpinBody
    s::Int
    E::Int
    i::Int
    j::Int
    k::Int

    SpinBody(i::Int, j::Int, N::Int) = new(rand([-1,1]), 0, i, j, k(i, j, N))
    SpinBody(s::Int, E::Int, i::Int, j::Int, k::Int) = new(s, E, i, j, k)
end

k(i, j, N) = i + (j - 1)*N

mutable struct SpinLattice
    n::Int
    N::Int
    bs::Matrix{SpinBody}
    T::Float64
    M::Float64
    E::Float64
    c::Float64
    steps::Int
    flips::Int
    f::Function

    function SpinLattice(n::Int, T::Float64; f::Function=amp)
        bs = Matrix{SpinBody}(undef, n, n)
        for i = 1:n, j = 1:n
            bs[i,j] = SpinBody(i, j, n)
        end
        init_energy!(bs, n)
        M = avg_magnetization(bs)
        E = avg_energy(bs)
        c = specific_heat(bs, T)
        new(n, n^2, bs, T, M, E, c, 0, 0, f)
    end

    SpinLattice(n::Int,
                N::Int,
                bs::Matrix{SpinBody},
                T::Float64,
                M::Float64,
                E::Float64,
                c::Float64,
                f::Function) = new(n, N, bs, T, M, E, c, 0, 0, f)
end

# TO-DO: add methods for correlation function

Base.copy(b::SpinBody) = SpinBody(b.s, b.E, b.i, b.j, b.k)

Base.copy(bs::Matrix{SpinBody}, n::Int) = begin
    bs_copy = Matrix{SpinBody}(undef, n, n);
    for k = eachindex(bs)
        bs_copy[k] = copy(bs[k])
    end;
    return bs_copy
end

Base.copy(L::SpinLattice) = SpinLattice(L.n, L.N, copy(L.bs, L.n), L.T, L.M, L.E, L.c, L.f)

avg_magnetization(bs::Matrix{SpinBody}) =
    mean([b.s for b in bs])

avg_energy(bs::Matrix{SpinBody}) =
    mean([b.E for b in bs])

specific_heat(bs::Matrix{SpinBody}, T::Float64) =
    (mean([b.E^2 for b in bs]) - mean([b.E for b in bs])^2) / T^2

function visualize(L::SpinLattice)
    vis = Matrix{Int}(undef, L.n, L.n)
    for k = eachindex(L.bs)
        vis[k] = L.bs[k].s
    end
    vis
end

function update_lattice!(L::SpinLattice)
    L.E = avg_energy(L.bs)
    L.M = avg_magnetization(L.bs)
    L.c = specific_heat(L.bs, L.T)
end

function init_energy!(bs::Matrix{SpinBody}, n::Int)
    for i = 1:n, j = 1:n
        sk = bs[i, j].s
        sr = bs[i, mod(j, n) + 1].s
        sl = bs[i, mod(j - 2, n) + 1].s
        su = bs[mod(i, n) + 1, j].s
        sd = bs[mod(i - 2, n) + 1, j].s
        bs[i,j].E = -sk * (sr + su + sl + sd)
    end
end

function flip!(L::SpinLattice, i::Int, j::Int; fast=true)
    L.bs[i,j].s *= -1
    n = L.n
    I = [mod(i, n) + 1, i, mod(i - 2, n) + 1, i]
    J = [j, mod(j, n) + 1, j, mod(j - 2, n) + 1]
    qs = collect(zip(I, J))
    if fast
        L.bs[i,j].E *= -1
        update_energies_fast!(L, qs, L.bs[i,j].s)
    else
        update_energies!(L, [(i, j); qs])
    end
end

function update_energies!(L::SpinLattice, qs::Vector{Tuple{Int,Int}})
    n = L.n
    for (i, j) in qs
        sk = L.bs[i, j].s
        sr = L.bs[i, mod(j, n) + 1].s
        sl = L.bs[i, mod(j - 2, n) + 1].s
        su = L.bs[mod(i, n) + 1, j].s
        sd = L.bs[mod(i - 2, n) + 1, j].s
        L.bs[i,j].E = -sk * (sr + su + sl + sd)
    end
end

function update_energies_fast!(L::SpinLattice, qs::Vector{Tuple{Int,Int}}, s::Int)
    for (i, j) in qs
        L.bs[i,j].s == s ? L.bs[i,j].E -= 2 : L.bs[i,j].E += 2
    end
end

prb(ΔE, T) = minimum([1, exp(-ΔE / T)])
amp(ΔE, T) = exp(-0.5 * ΔE / T)

function metropolis_step!(L::SpinLattice; fast=false)
    i, j = rand(1:L.N), rand(1:L.N)
    Ei = L.bs[i,j].E
    Ef = -Ei
    ΔE = Ef - Ei
    p = L.f(ΔE, L.T)
    u = rand()
    if u < p
        flip!(L, i, j, fast=fast)
        L.flips += 1
    end
    L.steps += 1
end

end
