# script to run the metropolis algorithm for the 2d ising model

using LinearAlgebra

mutable struct SpinBody
    s::Int64
    E::Float64
    i::Int64
    j::Int64
    k::Int64

    SpinBody(i::Int64, j::Int64) = new(rand([-1,1]), 0.0, i, j, (i - 1)*N + j)
end

# mutable struct PTree
#     sum::Float64
#     mid::Float64
#     l::Union{PTree, SpinBody, nothing}
#     r::Union{PTree, SpinBody, nothing}

#     PTree() = PTree(0.0, 0.0, nothing, nothing)
# end

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
            α = exp(-ΔE / T)
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
        println("C $(b.s == 1 ? 1 : 0) 0 $(b.s == 1 ? 0 : 1)")
        println("c $x $y $r")
    end
    println("F")
end

const T = 0.5
const N = 25

const steps = 1e5
const maxitr = 1e6

function main()
    bs = bodies()

    for step = 1:steps
        metropolis_step!(bs)
        anim(bs)
    end

    println("Q")
end

main()


