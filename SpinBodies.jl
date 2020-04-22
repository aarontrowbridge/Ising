#
# module for dynamic spin bodies
#

module SpinBodies

using Statistics

export SpinBody, SpinLattice, visualize, update_lattice!, update_energies!, metropolis_step!

mutable struct SpinBody
    s::Int64
    E::Float64
    i::Int64
    j::Int64
    k::Int64

    SpinBody(i::Int64, j::Int64, N::Int64) = new(rand([-1,1]), 0.0, i, j, k(i, j, N))
end

k(i, j, N) = (i - 1)*N + j

mutable struct SpinLattice
    N::Int64
    bs::Matrix{SpinBody}
    T::Float64
    M::Float64
    E::Float64
    c::Float64
    flips::Int64

    function SpinLattice(N::Int64, T::Float64)
        bs = Matrix{SpinBody}(undef, N, N)
        for i = 1:N, j = 1:N
            bs[i,j] = SpinBody(i, j, N)
        end
        init_energy!(bs, N)
        M = avg_magnetization(bs)
        E = avg_energy(bs)
        c = specific_heat(bs, T)
        new(N, bs, T, M, E, c, 0)
    end
end
# TO-DO: add methods for correlation function

avg_magnetization(bs::Matrix{SpinBody}) =
    mean([b.s for b in bs])

avg_energy(bs::Matrix{SpinBody}) =
    mean([b.E for b in bs])

specific_heat(bs::Matrix{SpinBody}, T::Float64) =
    (mean([b.E^2 for b in bs]) - mean([b.E for b in bs])^2) / T^2

function visualize(l::SpinLattice)
    vis = Matrix{Int64}(undef, l.N, l.N)
    for i = 1:N, j = 1:N
        vis[i,j] = l.bs[i,j].s
    end
    vis
end

function update_lattice!(l::SpinLattice)
    l.E = avg_energy(l.bs)
    l.M = avg_magnetization(l.bs)
    l.c = specific_heat(l.bs, l.T)
end

function init_energy!(bs::Matrix{SpinBody}, N::Int64)
    for i = 1:N, j = 1:N
        sk = bs[i, j].s
        sr = bs[i, mod(j, N) + 1].s
        su = bs[mod(i, N) + 1, j].s
        sl = bs[i, mod(j-2, N) + 1].s
        sd = bs[mod(i-2, N) + 1, j].s
        bs[i,j].E = -sk * (sr + su + sl + sd)
    end
end

function update_energies!(l::SpinLattice, i::Int64, j::Int64)
    l.bs[i,j].s *= -1
    l.bs[i,j].E *= -1
    N = l.N
    I = [mod(i, N) + 1, i, mod(i-2, N) + 1, i]
    J = [j, mod(i, N) + 1, j, mod(j-2, N) + 1]
    qs = collect(zip(I, J))
    update_energies!(l, qs)
end

function update_energies!(l::SpinLattice, qs::Vector{Tuple{Int64, Int64}})
    N = l.N
    for (i, j) in qs
        sk = l.bs[i, j].s
        sr = l.bs[i, mod(j, N) + 1].s
        su = l.bs[mod(i, N) + 1, j].s
        sl = l.bs[i, mod(j-2, N) + 1].s
        sd = l.bs[mod(i-2, N) + 1, j].s
        l.bs[i,j].E = -sk * (sr + su + sl + sd)
    end
end

m(ΔE, T) = minimum([1, exp(-ΔE / T)])

function metropolis_step!(l::SpinLattice)
    N = l.N
    i, j = rand(1:N), rand(1:N)
    Ei = l.bs[i,j].E
    Ef = -Ei
    ΔE = Ef - Ei
    p = m(ΔE, l.T)
    u = rand()
    if u < p
        l.flips += 1
        update_energies!(l, i, j)
    end
end

end
