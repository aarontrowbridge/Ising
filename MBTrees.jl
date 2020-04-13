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

    MBTree() = new(0.0, 0.0,
                   Vector{Int64}(undef, 0),
                   Vector{Int64}(undef, 0),
                   nothing, nothing, nothing, true)
end

function build_tree(bs, T)
    ptree = MBTree()
    for b in bs
        populate!(ptree, b, T)
    end
    ptree
end

m(ΔE, T) = minimum([1, exp(-ΔE / T)])

function populate!(t, b, T)
    ΔE = 2 * b.E
    p = m(ΔE, T)
    t.sum += p
    if t.b == nothing
        if t.leaf
            t.b = b
        else
            if rand(Bool)
                push!(t.lks, b.k)
                t.cut += p
                populate!(t.l, b, T)
            else
                push!(t.rks, b.k)
                populate!(t.r, b, T)
            end
        end
    else
        push!(t.lks, t.b.k)
        push!(t.rks, b.k)
        t.l = MBTree()
        t.r = MBTree()
        populate!(t.l, t.b, T)
        populate!(t.r, b, T)
        t.b = nothing
        t.leaf = false
    end
end

function choose_body(t)
    if t.leaf
        return (t.b.i, t.b.j)
    else
        x = rand(0:1e-10:t.sum)
        if x < t.cut
            choose_body(t.l)
        else
            choose_body(t.r)
        end
    end
end

function freeman_step!(bs, tree, N, T, maxitr)
    flip = false
    itr = 0
    while !flip
        (i, j) = choose_body(tree)

        Ei = bs[i,j].E
        Ef = -Ei

        if Ef < Ei
            bs[i,j].s *= -1
            tree_update_energies!(tree, bs, i, j, N, T)
            flip = true
        else
            ΔE = Ef - Ei
            α = m(ΔE, T)
            u = rand()
            if u < α
                bs[i,j].s *= -1
                tree_update_energies!(tree, bs, i, j, N, T)
                flip = true
            end
        end

        itr += 1
        if itr > maxitr
            println("Q")
            println("max iterations reached!")
            exit()
        end
    end
end

k(i, j, N) = (i - 1)*N + j

function tree_update_energies!(t, bs, i, j, N, T)
    I = [i, mod(i, N) + 1, i,   mod(i-2, N) + 1, i]
    J = [j, j,   mod(i, N) + 1, j,   mod(j-2, N) + 1]
    qs = collect(zip(I, J))
    update_energies!(bs, qs, N)

    ks = [k(i, j, N) for (i, j) in qs]
    ΔEs = energy_shifts(bs, qs, N)

    updates = collect(zip(ks, ΔEs))

    for (k, ΔE) in updates
        update!(t, k, ΔE, T)
    end
end

function energy_shifts(bs, qs, N)
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

function update!(t, k, ΔE, T)
    p = m(ΔE, T)
    if t.leaf
        t.sum += p
        t.b.E += ΔE
    else
        t.sum += p
        if k in t.lks
            t.cut += p
            update!(t.l, k, ΔE, T)
        else
            update!(t.r, k, ΔE, T)
        end
    end
end

end
