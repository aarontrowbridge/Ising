# script to run the metropolis algorithm for the 2d ising model

using LinearAlgebra

mutable struct SpinBody
    s::Int64
    E::Float64
    i::Int64
    j::Int64

    SpinBody(i::Int64, j::Int64) = new(rand([-1,1]), 0.0, i, j)
end

function bodies()
    bs = Matrix{SpinBody}(undef, N, N)

    for i = 1:N, j = 1:N
        bs[i,j] = SpinBody(i, j)
    end

    update_energy!(bs)

    bs
end

function update_energy!(bs)
    for i = 1:N, j = 1:N
        sk = bs[i, j].s

        sr = bs[i, mod(j, N) + 1].s
        su = bs[mod(i, N) + 1, j].s
        sl = bs[i, mod(j-2, N) + 1].s
        sd = bs[mod(i-2, N) + 1, j].s

        bs[i,j].E = -sk * (sr + su + sl + sd)
    end
end

function get_energy()

function metropolis_step!(bs)
    flipped = false

    while !flipped
        i, j = rand(1:N), rand(1:N)

        E = bs[i,j].E

        if -E < E
            bs[i,j].s *= -1
            update_energy(bs, i, j)
        else
            ΔE = 2E
            α = exp(-ΔE / T)
            u = rand()
            if u <= 


const T = 1
const N = 4
const steps = 1000

function main()
    bs = bodies()

    println(bs)

    # for step = 1:steps
    #     metropolis_step!(bs)
    #     anim(bs)
    # end

end

main()


