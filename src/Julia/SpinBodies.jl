#
# module for dynamic spin bodies
#

module SpinBodies

using Statistics
using WalterMethod

export SpinBody, SpinLattice
export ps, visualize, measure!, metropolis_step!, walter_step!

mutable struct SpinBody
    s::Int
    E::Int
    i::Int
    j::Int
    k::Int

    SpinBody(i::Int, j::Int, n::Int) = new(rand([-1,1]), 0, i, j, k(i, j, n))
    SpinBody(s::Int, E::Int, i::Int, j::Int, k::Int) = new(s, E, i, j, k)
end

k(i, j, n) = i + (j - 1) * n

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
    wsum::Float64
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
        new(n, n^2, bs, T, M, E, c, 0, 0, 0.0, f)
    end

    SpinLattice(n::Int,
                N::Int,
                bs::Matrix{SpinBody},
                T::Float64,
                M::Float64,
                E::Float64,
                c::Float64,
                f::Function) = new(n, N, bs, T, M, E, c, 0, 0, 0.0, f)
end

prb(ΔE, T) = minimum([1, exp(-ΔE / T)])
amp(ΔE, T) = exp(-0.5 * ΔE / T)

ps(L::SpinLattice) = vec([L.f(-2b.E, L.T) for b in L.bs])

# TO-DO: add methods for correlation function

Base.copy(b::SpinBody) = SpinBody(b.s, b.E, b.i, b.j, b.k)

Base.copy(bs::Matrix{SpinBody}, n::Int) = begin
    bs_copy = Matrix{SpinBody}(undef, n, n);
    for k = eachindex(bs)
        bs_copy[k] = copy(bs[k])
    end;
    return bs_copy
end

Base.copy(L::SpinLattice) = SpinLattice(
    L.n, L.N, copy(L.bs, L.n), L.T, L.M, L.E, L.c, L.f)

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
    return vis
end

function measure!(L::SpinLattice)
    L.E = avg_energy(L.bs)
    L.M = avg_magnetization(L.bs)
    L.c = specific_heat(L.bs, L.T)
end

@inline function init_energy!(bs::Matrix{SpinBody}, n::Int)
    for i = 1:n, j = 1:n
        sk = bs[i, j].s
        sr = bs[i, mod(j, n) + 1].s
        sl = bs[i, mod(j - 2, n) + 1].s
        su = bs[mod(i, n) + 1, j].s
        sd = bs[mod(i - 2, n) + 1, j].s
        bs[i,j].E = -sk * (sr + su + sl + sd)
    end
end

@inline function update_energies!(L::SpinLattice, qs::Vector{Tuple{Int,Int}}, s::Int)
    for (i, j) in qs
        L.bs[i,j].s == s ? L.bs[i,j].E -= 2 : L.bs[i,j].E += 2
    end
end

@inline function flip!(L::SpinLattice, i::Int, j::Int)
    n = L.n
    I = [mod(i, n) + 1, i, mod(i - 2, n) + 1, i]
    J = [j, mod(j, n) + 1, j, mod(j - 2, n) + 1]
    qs = collect(zip(I, J))
    L.bs[i,j].s *= -1
    L.bs[i,j].E *= -1
    update_energies!(L, qs, L.bs[i,j].s)
end

function metropolis_step!(L::SpinLattice)
    i, j = rand(1:L.n), rand(1:L.n)
    p = prb(-2*L.bs[i,j].E, L.T)
    u = rand()
    if u < p
        flip!(L, i, j)
        L.flips += 1
    end
    L.steps += 1
end

ij(k, n) = (mod(k - 1, n) + 1, div(k - 1, n) + 1)

function walter_step!(L::SpinLattice, tree::WalterTree)
    x = tree.psum * rand()
    n = L.n
    i, j = ij(move(tree, x), n)
    L.bs[i,j].s *= -1
    I = [mod(i, n) + 1, i, mod(i - 2, n) + 1, i]
    J = [j, mod(j, n) + 1, j, mod(j - 2, n) + 1]
    qs = collect(zip(I, J))
    Ei_s = [L.bs[x,y].E for (x, y) in [(i,j);qs]]
    L.bs[i,j].E *= -1
    update_energies!(L, qs, L.bs[i,j].s)
    Ef_s = [L.bs[x,y].E for (x, y) in [(i,j);qs]]
    ks = [k(x, y, n) for (x, y) in [(i,j);qs]]
    Δp_s = [L.f(-2 * Ef, L.T) - L.f(-2 * Ei, L.T) for (Ei, Ef) in zip(Ei_s, Ef_s)]
    update!(tree, collect(zip(ks, Δp_s))) # this causes code to break if n = 2
                                          # some ks are not unique => negatives in tree
    if L.f == prb
        Pₐ = tree.psum / L.N
        L.steps += Int(maximum([0, floor(log(1 - Pₐ, r))])) + 1
    else
        L.steps += 1
        L.wsum += L.N / tree.psum
    end
    L.flips += 1
end

end
