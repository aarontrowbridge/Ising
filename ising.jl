# script to run the metropolis algorithm for the 2d ising model

mutable struct SpinBody
    s::Int64
    E::Float64
    i::Int64
    j::Int64
    k::Int64

    SpinBody(i::Int64, j::Int64) = new(rand([-1,1]), 0.0, i, j, k(i, j))
end

k(i, j) = (i - 1)*N + j

mutable struct PTree
    sum::Float64
    cut::Float64
    lks::Vector{Int64}
    rks::Vector{Int64}
    b::Union{SpinBody, Nothing}
    l::Union{PTree, Nothing}
    r::Union{PTree, Nothing}
    leaf::Bool

    PTree() = new(0.0, 0.0,Vector{Int64}(undef, 0), Vector{Int64}(undef, 0), nothing, nothing, nothing, true)
end

function ptree(bs)
    ptree = PTree()
    for b in bs
        populate!(ptree, b)
    end
    ptree
end

p(ΔE) = minimum((1, exp(-ΔE / T)))

function populate!(t, b)
    ΔE = 2 * b.E
    pd = p(ΔE)
    t.sum += pd
    if t.b == nothing
        if t.leaf
            t.b = b
        else
            if rand(Bool)
                push!(t.lks, b.k)
                t.cut += pd
                populate!(t.l, b)
            else
                push!(t.rks, b.k)
                populate!(t.r, b)
            end
        end
    else
        push!(t.lks, t.b.k)
        push!(t.rks, b.k)
        t.l = PTree()
        t.r = PTree()
        populate!(t.l, t.b)
        populate!(t.r, b)
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



function freeman_step!(bs, tree)
    flipped = false

    itr = 0
    while !flipped
        (i, j) = choose_body(tree)

        Ei = bs[i,j].E
        Ef = -Ei

        if Ef < Ei
            bs[i,j].s *= -1
            tree_update_energies!(tree, bs, i, j)
            flipped = true
        else
            ΔE = Ef - Ei
            α = p(ΔE)
            u = rand()
            if u < α
                bs[i,j].s *= -1
                tree_update_energies!(tree, bs, i, j)
                flipped = true
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

function bodies()
    bs = Matrix{SpinBody}(undef, N, N)
    for i = 1:N, j = 1:N
        bs[i,j] = SpinBody(i, j)
    end
    init_energy!(bs)
    bs
end

function init_energy!(bs)
    for i = 1:N, j = 1:N
        sk = bs[i, j].s
        sr = bs[i, mod(j, N) + 1].s
        su = bs[mod(i, N) + 1, j].s
        sl = bs[i, mod(j-2, N) + 1].s
        sd = bs[mod(i-2, N) + 1, j].s
        bs[i,j].E = -sk * (sr + su + sl + sd)
    end
end

function update_adj_energies!(bs, i, j)
    I = [mod(i, N) + 1, i,   mod(i-2, N) + 1, i]
    J = [j,   mod(i, N) + 1, j,   mod(j-2, N) + 1]
    qs = collect(zip(I, J))
    update_energies!(bs, qs)
end

function tree_update_energies!(t, bs, i, j)
    I = [i, mod(i, N) + 1, i,   mod(i-2, N) + 1, i]
    J = [j, j,   mod(i, N) + 1, j,   mod(j-2, N) + 1]
    qs = collect(zip(I, J))
    update_energies!(bs, qs)

    ks = [k(i, j) for (i, j) in qs]
    ΔEs = energy_shifts(bs, qs)

    updates = collect(zip(ks, ΔEs))

    for (k, ΔE) in updates
        update!(t, k, ΔE)
    end
end

function update!(t, k, ΔE)
    pd = p(ΔE)
    if t.leaf
        t.sum += pd
        t.b.E += ΔE
    else
        t.sum += pd
        if k in t.lks
            t.cut += pd
            update!(t.l, k, ΔE)
        else
            update!(t.r, k, ΔE)
        end
    end
end

function energy_shifts(bs, qs)
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

function update_energies!(bs, qs)
    for (i, j) in qs
        sk = bs[i, j].s
        sr = bs[i, mod(j, N) + 1].s
        su = bs[mod(i, N) + 1, j].s
        sl = bs[i, mod(j-2, N) + 1].s
        sd = bs[mod(i-2, N) + 1, j].s
        bs[i,j].E = -sk * (sr + su + sl + sd)
    end
end

function metropolis_step!(bs)
    flipped = false

    itr = 0
    while !flipped
        i, j = rand(1:N), rand(1:N)

        Ei = bs[i,j].E
        Ef = -Ei

        if Ef < Ei
            bs[i,j].s *= -1
            bs[i,j].E = Ef
            update_adj_energies!(bs, i, j)
            flipped = true
        else
            ΔE = Ef - Ei
            α = p(ΔE)
            u = rand()
            if u < α
                bs[i,j].s *= -1
                bs[i,j].E = Ef
                update_adj_energies!(bs, i, j)
                flipped = true
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

function anim(bs)
    for b in bs
        x = b.i - N/2
        y = b.j - N/2
        r = 0.3
        println("C $(b.s == 1 ? 0 : 1) 0 $(b.s == 1 ? 1 : 0)")
        println("c $x $y $r")
    end
    println("F")
end

const T = 2.0
const N = 50

const steps = 1e8
const maxitr = 1e6

function main()
    bs = bodies()
    tree = ptree(bs)

    println(tree.sum)

    for step = 1:steps
        freeman_step!(bs, tree)
        # metropolis_step!(bs)
        if step % 50 == 0
            anim(bs)
        end
    end

    println("Q")
end

@time main()



