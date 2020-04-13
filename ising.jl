#
# script to run the metropolis algorithm for the 2d ising model
#

push!(LOAD_PATH, pwd())

using SpinBodies, MBTrees

function anim(bs)
    for b in bs
        x = b.i - N/2
        y = b.j - N/2
        r = 0.3
        println("C $(b.s == 1 ? 0 : 1) 0 $(b.s == 1 ? 1 : 0)")
        println("c3 $x $y 0 $r")
    end
    println("F")
end

const N = 40
const T = 2.0

const steps = 1e8
const maxitr = 1e6

function main()
    bs = spin_lattice(N)
    tree = build_tree(bs, T)

    for step = 1:steps

        freeman_step!(bs, tree, N, T, maxitr)
        # metropolis_step!(bs, N, T, maxitr)

        if step % 50 == 0 anim(bs) end
    end

    println("Q")
end

@time main()



