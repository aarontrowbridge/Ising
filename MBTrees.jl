#
# mudule for Metropolis Binary Trees
#

module MBTrees

using SpinBodies

export MBTree, build_tree, freeman_step!

mutable struct MBTree
    sum::Float64
    cut::Float64
    lks::Vector{Int64}
    rks::Vector{Int64}
    b::Union{SpinBody, Nothing}
    l::Union{MBTree, Nothing}
    r::Union{MBTree, Nothing}
    leaf::Bool
    flip::Bool
    N::Int64
    type::Symbol

    MBTree(N::Int64, type::Symbol) = new(0.0, 0.0,
                                         Vector{Int64}(undef, 0),
                                         Vector{Int64}(undef, 0),
                                         nothing, nothing, nothing,
                                         true, false,
                                         N, type)
end

function build_tree(l::SpinLattice, type::Symbol)
    N = l.N
    tree = MBTree(l.N, type)
    if type == :tri
        for b in l.bs
            populate_trinary!(tree, b, l.T)
        end
    elseif type == :bi
        for b in l.bs
            populate_binary!(tree, b, l.T)
        end
    end
    tree
end

m(ΔE, T) = minimum([1, exp(-ΔE / T)])

function populate_trinary!(t::MBTree, b::SpinBody, T::Float64)
    p = m(2*b.E, T)
    t.sum += m(2*b.E, T)
    if t.b == nothing
        t.b = b
    else
        t.leaf = false
        if rand(Bool)
            push!(t.lks, b.k)
            if t.l == nothing
                t.l = MBTree(t.N, t.type)
            end
            populate_trinary!(t.l, b, T)
        else
            push!(t.rks, b.k)
            if t.r == nothing
                t.r = MBTree(t.N, t.type)
            end
            populate_trinary!(t.r, b, T)
        end
    end
end

function populate_binary!(t::MBTree, b::SpinBody, T::Float64)
    p = m(2*b.E, T)
    t.sum += p
    if t.b == nothing
        if t.leaf
            t.b = b
        else
            if rand(Bool)
                push!(t.lks, b.k)
                t.cut += p
                populate_binary!(t.l, b, T)
            else
                push!(t.rks, b.k)
                populate_binary!(t.r, b, T)
            end
        end
    else
        push!(t.lks, t.b.k)
        push!(t.rks, b.k)
        t.l = MBTree(t.N, t.type)
        t.r = MBTree(t.N, t.type)
        populate_binary!(t.l, t.b, T)
        populate_binary!(t.r, b, T)
        t.b = nothing
        t.leaf = false
    end
end

function choose_body_trinary(t::MBTree, x::Float64)
    if t.leaf
        return (t.b.i, t.b.j)
    else
        if t.r == nothing
            if x < t.l.sum
                choose_body_trinary(t.l, x)
            else
                return (t.b.i, t.b.j)
            end
        elseif t.l == nothing
            if x > t.sum - t.r.sum
                choose_body_trinary(t.r, t.sum - x)
            else
                return (t.b.i, t.b.j)
            end
        else
            cut1 = t.l.sum
            cut2 = t.sum - t.r.sum
            if x < cut1
                choose_body_trinary(t.l, x)
            elseif x > cut2
                choose_body_trinary(t.r, t.sum - x)
            else
                return (t.b.i, t.b.j)
            end
        end
    end
end

function choose_body_binary(t::MBTree, x::Float64)
    if t.leaf
        return (t.b.i, t.b.j)
    else
        if x < t.cut
            choose_body_binary(t.l, x)
        else
            choose_body_binary(t.r, t.sum - x)
        end
    end
end

function freeman_step!(l::SpinLattice, tree::MBTree)
    x = rand(0:1e-10:tree.sum)
    if tree.type == :tri
        (i, j) = choose_body_trinary(tree, x)
    else
        (i, j) = choose_body_binary(tree, x)
    end
    Ei = l.bs[i,j].E
    Ef = -Ei
    ΔE = Ef - Ei
    α = m(ΔE, l.T)
    u = rand()
    if u < α
        l.bs[i,j].s *= -1
        tree_update_energies!(tree, l, i, j)
        tree.flip = true
        l.flips += 1
    end
end

k(i, j, N) = (i - 1)*N + j

function tree_update_energies!(tree::MBTree, l::SpinLattice, i::Int64, j::Int64)
    N = l.N
    l.bs[i,j].E *= -1
    I = [mod(i, N) + 1, i,   mod(i-2, N) + 1, i]
    J = [j,   mod(i, N) + 1, j,   mod(j-2, N) + 1]
    qs = collect(zip(I, J))

    update_energies!(l, qs)

    ks = [k(i, j, N) for (i, j) in qs]
    ΔEs = energy_shifts(l.bs, qs)

    if tree.type == :tri
        for (k, ΔE) in zip(ks, ΔEs)
            p = m(ΔE, l.T)
            update_trinary!(tree, k, p, ΔE, l.T)
        end
    else
        for (k, ΔE) in zip(ks, ΔEs)
            p = m(ΔE, l.T)
            update_binary!(tree, k, p, ΔE, l.T)
        end
    end
end

function energy_shifts(bs::Matrix{SpinBody}, qs::Vector{Tuple{Int64,Int64}})
    N = size(bs)[1]
    ΔEs = []
    for (i, j) in qs
        sk = bs[i, j].s
        sr = bs[i, mod(j, N) + 1].s
        su = bs[mod(i, N) + 1, j].s
        sl = bs[i, mod(j-2, N) + 1].s
        sd = bs[mod(i-2, N) + 1, j].s
        Ef = -sk * (sr + su + sl + sd)
        Ei = bs[i,j].E
        ΔE = Ef - Ei
        push!(ΔEs, ΔE)
    end
    ΔEs
end

function update_trinary!(t::MBTree, k::Int64, p::Float64, ΔE::Float64, T::Float64)
    t.sum += p
    if t.b.k == k
        t.b.E += ΔE
    else
        if k in t.lks
            update_trinary!(t.l, k, p, ΔE, T)
        else
            update_trinary!(t.r, k, p, ΔE, T)
        end
    end
end


function update_binary!(t::MBTree, k::Int64, p::Float64, ΔE::Float64, T::Float64)
    t.sum += p
    if t.leaf
        t.b.E += ΔE
    else
        if k in t.lks
            t.cut += p
            update_binary!(t.l, k, p, ΔE, T)
        else
            update_binary!(t.r, k, p, ΔE, T)
        end
    end
end

end
