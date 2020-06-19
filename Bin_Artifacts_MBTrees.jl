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
    lks::Vector{Int64}
    rks::Vector{Int64}
    l::Union{MBTree, Nothing}
    r::Union{MBTree, Nothing}
    leaf::Bool
    N::Int64
    type::Symbol

    MBTree(N::Int64, type::Symbol) = new((i=0, j=0, k=0), 0., 0.,
                                         Vector{Int64}(undef, 0),
                                         Vector{Int64}(undef, 0),
                                         nothing, nothing, true,
                                         N, type)
end

# macro populate(type::Symbol)
#     eval(populate_$(string(type))nary!(tree, b))

function build_tree(l::SpinLattice, type::Symbol)
    N = l.N
    tree = MBTree(l.N, type)
    if type == :tri
        for b in l.bs
            p = m(-2*b.E, l.T)
            populate_trinary!(tree, b, p)
        end
    elseif type == :bi
        for b in l.bs
            p = m(-2*b.E, l.T)
            populate_binary!(tree, b, p)
        end
    end
    tree
end

m(ΔE, T) = minimum([1, exp(-ΔE / T)])

function populate_trinary!(t::MBTree, b::SpinBody, p::Float64)
    t.sum += p
    if t.p == 0.
        t.p = p
        t.q = (i=b.i, j=b.j, k=b.k)
    else
        t.leaf = false
        if rand(Bool)
            push!(t.lks, b.k)
            if t.l == nothing
                t.l = MBTree(t.N, t.type)
            end
            populate_trinary!(t.l, b, p)
        else
            push!(t.rks, b.k)
            if t.r == nothing
                t.r = MBTree(t.N, t.type)
            end
            populate_trinary!(t.r, b, p)
        end
    end
end

function populate_binary!(t::MBTree, b::SpinBody, p::Float64)
    t.sum += p
    if t.p == 0.
        if t.leaf
            t.p = p
            t.q = (i=b.i, j=b.j, k=b.k)
        else
            if rand(Bool)
                push!(t.lks, b.k)
                populate_binary!(t.l, b, p)
            else
                push!(t.rks, b.k)
                populate_binary!(t.r, b, p)
            end
        end
    else
        push!(t.lks, t.q.k)
        push!(t.rks, b.k)
        t.l = MBTree(t.N, t.type)
        t.r = MBTree(t.N, t.type)
        t.l.p = t.p
        t.p = 0.
        t.l.q = t.q
        t.q = (i=0, j=0, k=0)
        populate_binary!(t.r, b, p)
        t.leaf = false
    end
end

function choose_body_trinary(t::MBTree, x::Float64)
    if t.leaf
        return t.q
    else
        if t.r == nothing
            if x < t.l.sum
                choose_body_trinary(t.l, x)
            else
                return t.q
            end
        elseif t.l == nothing
            if x > t.p
                choose_body_trinary(t.r, x - t.p)
            else
                return t.q
            end
        else
            cut1 = t.l.sum
            cut2 = cut1 + t.p
            if x < cut1
                choose_body_trinary(t.l, x)
            elseif x > cut2
                choose_body_trinary(t.r, x - cut2)
            else
                return t.q
            end
        end
    end
end

function choose_body_binary(t::MBTree, x::Float64)
    if t.leaf
        return t.q
    else
        if x < t.l.sum
            choose_body_binary(t.l, x)
        else
            choose_body_binary(t.r, x - t.l.sum)
        end
    end
end

function freeman_step!(l::SpinLattice, t::MBTree)
    x = rand() * t.sum
    q = choose_body_trinary(t, x)
    tree_flip!(t, l, q.i, q.j)
end

function freeman_step_fast!(l::SpinLattice, t::MBTree)
    x = rand() * t.sum
    q = choose_body_trinary(t, x)
    tree_flip_fast!(t, l, q.i, q.j)
end

k(i, j, N) = (i - 1)*N + j

function tree_flip!(t::MBTree, l::SpinLattice, i::Int64, j::Int64)
    N = l.N
    I = [mod(i, N) + 1, i, mod(i-2, N) + 1, i]
    J = [j, mod(i, N) + 1, j, mod(j-2, N) + 1]
    qs = collect(zip(I, J))
    Ei_s = [l.bs[n,m].E for (n, m) in [(i,j);qs]]
    l.bs[i,j].s *= -1
    l.bs[i,j].E *= -1
    update_energies!(l, qs)
    Ef_s = [l.bs[n,m].E for (n, m) in [(i,j);qs]]
    ks = [k(n, m, N) for (n, m) in [(i,j);qs]]
    Δp_s = [m(-2*Ef, l.T) - m(-2*Ei, l.T) for (Ei, Ef) in zip(Ei_s, Ef_s)]
    for (k, Δp) in zip(ks, Δp_s)
        update_ps!(t, k, Δp)
    end
    l.flips += 1
end

function tree_flip_fast!(t::MBTree, l::SpinLattice, i::Int64, j::Int64)
    N = l.N
    I = [mod(i, N) + 1, i, mod(i-2, N) + 1, i]
    J = [j, mod(i, N) + 1, j, mod(j-2, N) + 1]
    qs = collect(zip(I, J))
    Ei_s = [l.bs[n,m].E for (n, m) in [(i,j);qs]]
    l.bs[i,j].s *= -1
    l.bs[i,j].E *= -1
    update_energies_fast!(l, qs, l.bs[i,j].s)
    Ef_s = [l.bs[n,m].E for (n, m) in [(i,j);qs]]
    ks = [k(n, m, N) for (n, m) in [(i,j);qs]]
    Δp_s = [m(-2*Ef, l.T) - m(-2*Ei, l.T) for (Ei, Ef) in zip(Ei_s, Ef_s)]
    for (k, Δp) in zip(ks, Δp_s)
        update_ps!(t, k, Δp)
    end
    l.flips += 1
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

function update_ps!(t::MBTree, k::Int64, Δp::Float64)
    t.sum += Δp
    if t.q.k == k
        t.p += Δp
    else
        if k in t.lks
            update_ps!(t.l, k, Δp)
        else
            update_ps!(t.r, k, Δp)
        end
    end
end

end
