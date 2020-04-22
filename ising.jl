#
# script to run the metropolis algorithm for the 2d ising model
#

using Printf

push!(LOAD_PATH, pwd())

using SpinBodies, MBTrees

function anim(bs, step, flip)
    for b in bs
        x = b.i - N/2
        y = b.j - N/2
        r = 0.4
        @printf "C %d 0 %d\n" (b.s == 1 ? 0 : 1) (b.s == 1 ? 1 : 0)
        @printf "c3 %f %f 0 %f\n" x y r
    end

    @printf "C 1 1 1\n"
    @printf "T -0.95 0.95\n"
    @printf "step = %d\n" step
    @printf "T -0.35 0.95\n"
    @printf "flip = %d\n" flip
    @printf "T 0.35 0.95\n"
    @printf "f/s = %.4f\n" flip/step
    @printf "F\n"
end

const N = 50
const T = 0.5

const steps = 1e7
const maxitr = 1e6

const type = :tri

const freeman = true

function main()
    bs = spin_lattice(N)

    if freeman
        tree = build_tree(bs, T, type)
    end

    flips = 0
    itr = 0

    for step = 1:steps
        if freeman
            freeman_step!(bs, tree)
            if tree.flip
                itr = 0
                flips += 1
                tree.flip = false
            end
        else
            (bs, flip) = metropolis_step!(bs, T)
            if flip
                itr = 0
                flips += 1
            end
        end

        if step % 500 == 0
            anim(bs, step, flips)
        end

        if itr > maxitr
            println("Q")
            println("reached max iterations")
            exit()
        end
        itr += 1
    end

    println("Q")
end

@time main()



