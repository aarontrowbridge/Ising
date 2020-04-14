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
    T::Float64
    N::Int64
    type::Symbol

    MBTree(T::Float64, N::Int64, type::Symbol) = new(0.0, 0.0,
                                                     Vector{Int64}(undef, 0),
                                                     Vector{Int64}(undef, 0),
                                                     nothing, nothing, nothing,
                                                     true, false,
                                                     T, N, type)
end

function build_tree(bs, T, type)
    N = size(bs)[1]
    tree = MBTree(T, N, type)
    if type == :tri
        for b in bs
            populate_trinary!(tree, b)
        end
    elseif type == :bi
        for b in bs
            populate_binary!(tree, b)
        end
    end
    tree
end

m(ΔE, T) = minimum([1, exp(-ΔE / T)])

function populate_trinary!(t, b)
    t.sum += m(2*b.E, t.T)
    if t.b == nothing
        t.b = b
    else
        t.leaf = false
        if rand(Bool)
            push!(t.lks, b.k)
            if t.l == nothing
                t.l = MBTree(t.T, t.N, t.type)
            end
            populate_trinary!(t.l, b)
        else
            push!(t.rks, b.k)
            if t.r == nothing
                t.r = MBTree(t.T, t.N, t.type)
            end
            populate_trinary!(t.r, b)
        end
    end
end

function choose_body_trinary(t, x)
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

function populate_binary!(t, b)
    p = m(2*b.E, t.T)
    t.sum += p
    if t.b == nothing
        if t.leaf
            t.b = b
        else
            if rand(Bool)
                push!(t.lks, b.k)
                t.cut += p
                populate_binary!(t.l, b)
            else
                push!(t.rks, b.k)
                populate_binary!(t.r, b)
            end
        end
    else
        push!(t.lks, t.b.k)
        push!(t.rks, b.k)
        t.l = MBTree(t.T, t.N, t.type)
        t.r = MBTree(t.T, t.N, t.type)
        populate_binary!(t.l, t.b)
        populate_binary!(t.r, b)
        t.b = nothing
        t.leaf = false
    end
end

function choose_body_binary(t, x)
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

function freeman_step!(bs, tree)
    x = rand(0:1e-10:tree.sum)
    if tree.type == :tri
        (i, j) = choose_body_trinary(tree,x)
    else
        (i, j) = choose_body_binary(tree, x)
    end
    Ei = bs[i,j].E
    Ef = -Ei
    ΔE = Ef - Ei
    α = m(ΔE, tree.T)
    u = rand()
    if u < α
        bs[i,j].s *= -1
        tree_update_energies!(tree, bs, i, j)
        tree.flip = true
    end
end

k(i, j, N) = (i - 1)*N + j

function tree_update_energies!(tree, bs, i, j)
    N = tree.N
    I = [i, mod(i, N) + 1, i,   mod(i-2, N) + 1, i]
    J = [j, j,   mod(i, N) + 1, j,   mod(j-2, N) + 1]
    qs = collect(zip(I, J))
    update_energies!(bs, qs)

    ks = [k(i, j, N) for (i, j) in qs]
    ΔEs = energy_shifts(bs, qs)
    if tree.type == :tri
        for (k, ΔE) in zip(ks, ΔEs)
            update_trinary!(tree, k, ΔE)
        end
    else
        for (k, ΔE) in zip(ks, ΔEs)
            update_binary!(tree, k, ΔE)
        end
    end
end

function energy_shifts(bs, qs)
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

function update_trinary!(t, k, ΔE)
    p = m(ΔE, t.T)
    t.sum += p
    if t.b.k == k
        t.b.E += ΔE
    else
        if k in t.lks
            update_trinary!(t.l, k, ΔE)
        else
            update_trinary!(t.r, k, ΔE)
        end
    end
end

function update_binary!(t, k, ΔE)
    p = m(ΔE, t.T)
    t.sum += p
    if t.leaf
        t.b.E += ΔE
    else
        if k in t.lks
            t.cut += p
            update_binary!(t.l, k, ΔE)
        else
            update_binary!(t.r, k, ΔE)
        end
    end
end

end
