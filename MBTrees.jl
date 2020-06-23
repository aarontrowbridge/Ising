#
# mudule for Metropolis Binary Trees
#

module MBTrees

using SpinBodies

export MBTree, build_tree, freeman_step!, freeman_step_fast!

mutable struct MBTree
    q::NamedTuple{(:i, :j, :k),Tuple{Int64, Int64, Int64}}
    p::Float64
    sum::Float64
    left_ks::Vector{Int64}
    right_ks::Vector{Int64}
    left::Union{MBTree, Nothing}
    right::Union{MBTree, Nothing}
    leaf::Bool
    N::Int64

    MBTree(N::Int64) = new((i=0, j=0, k=0), 0., 0.,
                           Vector{Int64}(undef, 0),
                           Vector{Int64}(undef, 0),
                           nothing, nothing, true, N)
end

function build_tree(L::SpinLattice)
    N = L.N
    tree = MBTree(L.N)
    for b in L.bs
        p = m(-2*b.E, L.T)
        populate!(tree, b, p)
    end
    tree
end

m(ΔE, T) = minimum([1, exp(-ΔE / T)])

function populate!(T::MBTree, b::SpinBody, p::Float64)
    T.sum += p
    if T.p == 0.
        T.p = p
        T.q = (i=b.i, j=b.j, k=b.k)
    else
        T.leaf = false
        if rand(Bool)
            push!(T.left_ks, b.k)
            if T.left == nothing
                T.left = MBTree(T.N)
            end
            populate!(T.left, b, p)
        else
            push!(T.right_ks, b.k)
            if T.right == nothing
                T.right = MBTree(T.N)
            end
            populate!(T.right, b, p)
        end
    end
end

function choose_body(T::MBTree, x::Float64)
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

function freeman_step!(L::SpinLattice, T::MBTree)
    x = rand() * T.sum
    b = choose_body(T, x)
    tree_flip!(L, T, b.i, b.j)
end

function freeman_step_fast!(L::SpinLattice, T::MBTree)
    x = rand() * T.sum
    b = choose_body(T, x)
    tree_flip_fast!(L, T, b.i, b.j)
end

k(i, j, N) = (i - 1)*N + j

function tree_flip!(L::SpinLattice, T::MBTree, i::Int64, j::Int64)
    N = L.N
    I = [mod(i, N) + 1, i, mod(i-2, N) + 1, i]
    J = [j, mod(j, N) + 1, j, mod(j-2, N) + 1]
    qs = collect(zip(I, J))
    Ei_s = [L.bs[n,m].E for (n, m) in [(i,j);qs]]
    L.bs[i,j].s *= -1
    L.bs[i,j].E *= -1
    update_energies!(L, qs)
    Ef_s = [L.bs[n,m].E for (n, m) in [(i,j);qs]]
    ks = [k(n, m, N) for (n, m) in [(i,j);qs]]
    Δp_s = [m(-2 * Ef, L.T) - m(-2 * Ei, L.T) for (Ei, Ef) in zip(Ei_s, Ef_s)]
    for (k, Δp) in zip(ks, Δp_s)
        update_ps!(T, k, Δp)
    end
    L.flips += 1
end

function update_energies!(L::SpinLattice, qs::Vector{Tuple{Int,Int}})
    N = L.N
    for (i, j) in qs
        sk = L.bs[i, j].s
        sr = L.bs[i, mod(j, N) + 1].s
        sl = L.bs[i, mod(j-2, N) + 1].s
        su = L.bs[mod(i, N) + 1, j].s
        sd = L.bs[mod(i-2, N) + 1, j].s
        L.bs[i,j].E = -sk * (sr + sl + su + sd)
    end
end

function tree_flip_fast!(L::SpinLattice, T::MBTree, i::Int64, j::Int64)
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
    Δp_s = [m(-2 * Ef, L.T) - m(-2 * Ei, L.T) for (Ei, Ef) in zip(Ei_s, Ef_s)]
    for (k, Δp) in zip(ks, Δp_s)
        update_ps!(T, k, Δp)
    end
    L.flips += 1
end

function update_energies_fast!(L::SpinLattice, qs::Vector{Tuple{Int,Int}}, s::Int)
    for (i, j) in qs
        if L.bs[i,j].s == s
            L.bs[i,j].E -= 2
        else
            L.bs[i,j].E += 2
        end
    end
end

function update_ps!(T::MBTree, k::Int64, Δp::Float64)
    T.sum += Δp
    if T.q.k == k
        T.p += Δp
    else
        if k in T.left_ks
            update_ps!(T.left, k, Δp)
        else
            update_ps!(T.right, k, Δp)
        end
    end
end

end
