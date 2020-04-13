#
# module for dynamic spin bodies
#

module SpinBodies

export SpinBody, spin_lattice, update_energies!, metropolis_step!

mutable struct SpinBody
    s::Int64
    E::Float64
    i::Int64
    j::Int64
    k::Int64

    SpinBody(i::Int64, j::Int64, N::Int64) = new(rand([-1,1]), 0.0,
                                                 i, j, k(i, j, N))
end

k(i, j, N) = (i - 1)*N + j

function init_energy!(bs, N)
    for i = 1:N, j = 1:N
        sk = bs[i, j].s
        sr = bs[i, mod(j, N) + 1].s
        su = bs[mod(i, N) + 1, j].s
        sl = bs[i, mod(j-2, N) + 1].s
        sd = bs[mod(i-2, N) + 1, j].s
        bs[i,j].E = -sk * (sr + su + sl + sd)
    end
end

function spin_lattice(N)
    bs = Matrix{SpinBody}(undef, N, N)
    for i = 1:N, j = 1:N
        bs[i,j] = SpinBody(i, j, N)
    end
    init_energy!(bs, N)
    bs
end

function update_adj_energies!(bs, i, j, N)
    I = [mod(i, N) + 1, i,   mod(i-2, N) + 1, i]
    J = [j,   mod(i, N) + 1, j,   mod(j-2, N) + 1]
    qs = collect(zip(I, J))
    update_energies!(bs, qs, N)
end

function update_energies!(bs, qs, N)
    for (i, j) in qs
        sk = bs[i, j].s
        sr = bs[i, mod(j, N) + 1].s
        su = bs[mod(i, N) + 1, j].s
        sl = bs[i, mod(j-2, N) + 1].s
        sd = bs[mod(i-2, N) + 1, j].s
        bs[i,j].E = -sk * (sr + su + sl + sd)
    end
end

m(ΔE, T) = minimum((1, exp(-ΔE / T)))

function metropolis_step!(bs, N, T, maxitr)
    flip = false
    itr = 0
    while !flip
        i, j = rand(1:N), rand(1:N)

        Ei = bs[i,j].E
        Ef = -Ei

        if Ef < Ei
            bs[i,j].s *= -1
            bs[i,j].E = Ef
            update_adj_energies!(bs, i, j, N)
            flip = true
        else
            ΔE = Ef - Ei
            α = m(ΔE, T)
            u = rand()
            if u < α
                bs[i,j].s *= -1
                bs[i,j].E = Ef
                update_adj_energies!(bs, i, j, N)
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


end
