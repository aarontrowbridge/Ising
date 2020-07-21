#
# mudule for Metropolis Binary Trees
#

module MBTrees

using SpinBodies

export MBTree, TreeVec
export build_tree, freeman_step!, treevecstep!

const Flt = Float64

mutable struct TreeVec
    cps::Vector{Flt}
    sum::Flt

    function TreeVec(L::SpinLattice)
        ps = [L.f(-2 * b.E, L.T) for b in L.bs]
        cps = [sum(ps[1:k]) for k = 1:(L.N^2 - 1)]
        sm = sum(ps)
        new(cps, sm)
    end
end

function treevecstep!(L::SpinLattice, treevec::TreeVec)
    r = treevec.sum * rand()
    k = searchsorted(treevec.cps, r).start
    i, j = ij(k, L.N)
    updates = treevecflip!(L, i, j)
    update!(treevec, updates, L.N)
    if L.f == prb
        Pₐ = treevec.sum / L.N^2
        L.steps += Int(maximum([0, floor(log(1 - Pₐ, r))])) + 1
    else
        L.steps += 1
    end
    L.flips += 1
end

ij(k, N) = (div(k - 1, N) + 1, mod(k - 1, N) + 1)

function update!(treevec::TreeVec, updates::Vector, N::Int)
    for (k, Δp) in updates
        if k < N^2 treevec.cps[k:end] .+= Δp end
        treevec.sum += Δp
    end
end

function treevecflip!(L::SpinLattice, i::Int, j::Int)
    N = L.N
    I = [mod(i, N) + 1, i, mod(i-2, N) + 1, i]
    J = [j, mod(j, N) + 1, j, mod(j-2, N) + 1]
    qs = collect(zip(I, J))
    Ei_s = [L.bs[n,m].E for (n, m) in [(i,j);qs]]
    L.bs[i,j].s *= -1
    L.bs[i,j].E *= -1
    update_energies_fast!(L, qs, L.bs[i,j].s)
    Ef_s = [L.bs[n,m].E for (n, m) in [(i,j);qs]]
    ks = [k(n, m, N) for (n, m) in [(i,j);qs]]
    Δp_s = [L.f(-2 * Ef, L.T) - L.f(-2 * Ei, L.T) for (Ei, Ef) in zip(Ei_s, Ef_s)]
    return collect(zip(ks, Δp_s))
end

mutable struct MBTree
    q::NamedTuple{(:i,:j,:k),Tuple{Int,Int,Int}}
    p::Union{Flt,Nothing}
    sum::Flt
    left_ks::Vector{Int}
    right_ks::Vector{Int}
    up::Union{MBTree, Nothing}
    left::Union{MBTree,Nothing}
    right::Union{MBTree,Nothing}
    leaf::Bool

    MBTree() = new((i=0, j=0, k=0), nothing, 0.,
                   Vector{Int}(undef, 0),
                   Vector{Int}(undef, 0),
                   nothing, nothing, nothing, true)

    MBTree(tree::MBTree) = new((i=0, j=0, k=0), nothing, 0.,
                               Vector{Int}(undef, 0),
                               Vector{Int}(undef, 0),
                               tree, nothing, nothing, true)

end

function build_tree(L::SpinLattice)
    tree = MBTree()
    for b in L.bs
        p = L.f(-2*b.E, L.T)
        populate!(tree, b, Flt(p))
    end
    tree
end

function populate!(T::MBTree, b::SpinBody, p::Flt)
    T.sum += p
    if T.p == nothing
        T.p = p
        T.q = (i=b.i, j=b.j, k=b.k)
    else
        T.leaf = false
        if rand(Bool)
            push!(T.left_ks, b.k)
            if T.left == nothing
                T.left = MBTree(T)
            end
            populate!(T.left, b, p)
        else
            push!(T.right_ks, b.k)
            if T.right == nothing
                T.right = MBTree(T)
            end
            populate!(T.right, b, p)
        end
    end
end

function choose_body(T::MBTree, x)
    if T.leaf
        return T.q
    else
        if T.right == nothing
            if x < T.left.sum
                choose_body(T.left, x)
            else
                return T.q
            end
        elseif T.left == nothing
            if x > T.p
                choose_body(T.right, x - T.p)
            else
                return T.q
            end
        else
            cut1 = T.left.sum
            cut2 = cut1 + T.p
            if x < cut1
                choose_body(T.left, x)
            elseif x > cut2
                choose_body(T.right, x - cut2)
            else
                return T.q
            end
        end
    end
end

function freeman_step!(L::SpinLattice, T::MBTree; fast=false)
    r = rand(Flt)
    b = choose_body(T, r * T.sum)
    tree_flip!(L, T, b.i, b.j, r, fast=fast)
end

k(i, j, N) = (i - 1)*N + j

function tree_flip!(L::SpinLattice, T::MBTree, i::Int, j::Int, r::Flt; fast=false)
    N = L.N
    I = [mod(i, N) + 1, i, mod(i-2, N) + 1, i]
    J = [j, mod(j, N) + 1, j, mod(j-2, N) + 1]
    qs = collect(zip(I, J))
    Ei_s = [L.bs[n,m].E for (n, m) in [(i,j);qs]]
    L.bs[i,j].s *= -1
    if fast
        L.bs[i,j].E *= -1
        update_energies_fast!(L, qs, L.bs[i,j].s)
    else
        update_energies!(L, [(i,j);qs])
    end
    Ef_s = [L.bs[n,m].E for (n, m) in [(i,j);qs]]
    ks = [k(n, m, N) for (n, m) in [(i,j);qs]]
    Δp_s = [L.f(-2 * Ef, L.T) - L.f(-2 * Ei, L.T) for (Ei, Ef) in zip(Ei_s, Ef_s)]
    for (k, Δp) in zip(ks, Δp_s)
        update_tree!(T, k, Flt(Δp))
    end
    L.flips += 1
    if L.f == prb
        Pₐ = T.sum / N^2
        L.steps += Int(maximum([0, floor(log(1 - Pₐ, r))])) + 1
    else
        L.steps += 1
    end
end

function update_tree!(T::MBTree, k::Int, Δp::Flt)
    T.sum += Δp
    if T.q.k == k
        T.p += Δp
    else
        if k in T.left_ks
            update_tree!(T.left, k, Δp)
        else
            update_tree!(T.right, k, Δp)
        end
    end
end




end
