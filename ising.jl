#
# script to run the metropolis algorithm for the 2d ising model
#

using Printf

push!(LOAD_PATH, pwd())

using SpinBodies, MBTrees

function anim(l, steps, flips)
    for b in l.bs
        x = b.i - N/2
        y = b.j - N/2
        r = 0.4
        @printf "C %d 0 %d\n" (b.s == 1 ? 0 : 1) (b.s == 1 ? 1 : 0)
        @printf "c3 %f %f 0 %f\n" x y r
    end

    @printf "C 1 1 1\n"
    @printf "T -0.95 0.95\n"
    @printf "steps = %d\n" steps
    @printf "T -0.35 0.95\n"
    @printf "flips = %d\n" l.flips
    @printf "T 0.25 0.95\n"
    @printf "steps/flip = %.4f\n" frameskip/maximum([1, l.flips - flips])
    @printf "F\n"
end


const steps = 1e6
const maxitr = 1e6

const frameskip = 500

const freeman = false
const type = :bi

const N = 50
const T = 0.5

function main()
    lattice = SpinLattice(N, T)

    if freeman
        tree = build_tree(lattice, type)
    end

    flips = 0
    for step = 1:steps
        if freeman
            freeman_step!(lattice, tree)
        else
            metropolis_step!(lattice)
        end

        if step % frameskip == 0
            anim(lattice, step, flips)
            flips = lattice.flips
        end
    end

    println("Q")
end

@time main()



